import os
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
import numpy as np
import moments

# -----------------------
# USER SETTINGS
# -----------------------
SYNC_DIR = "."          # folder containing *.masked.sync.gz
OUT_DIR = "per_gene_SFS_from_sync_mindepth"
GENE_BED = "fbgn_genes.bed"

FOURFOLD_BED = "fourfold_sites.destnames.bed"
ZEROFOLD_BED = "zerofold_sites.destnames.bed"

KEEP_CHROMS = {"2L", "2R", "3L", "3R", "X", "4", "Y"}

n_proj = 20
ALLOW_MULTI_ALLELIC = False  #skip sites with >1 non-ref allele >0

# minimum called depth (A+T+C+G) to consider a site callable
# (callable if A+T+C+G => 10). This choice is from https://academic.oup.com/mbe/article/38/12/5782/6361628
MIN_DEPTH = 10

# -----------------------
# HELPERS
# -----------------------
BASE_TO_IDX = {"A": 0, "T": 1, "C": 2, "G": 3}


def load_genes_bed(gene_bed: str) -> Dict[str, List[Tuple[int, int, str]]]:
    genes: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    with open(gene_bed, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, gid = line.rstrip("\n").split("\t")[:4]
            if chrom not in KEEP_CHROMS:
                continue
            genes[chrom].append((int(start), int(end), gid))
    for chrom in genes:
        genes[chrom].sort(key=lambda x: x[0])
    return genes


def find_gene_for_pos0(
    genes_chr: List[Tuple[int, int, str]],
    pos0: int,
    ptr: int,
) -> Tuple[Optional[str], int]:
    n = len(genes_chr)
    while ptr < n and genes_chr[ptr][1] <= pos0:
        ptr += 1

    i = ptr
    while i < n and genes_chr[i][0] <= pos0:
        start0, end1, gid = genes_chr[i]
        if pos0 < end1:
            return gid, ptr
        i += 1

    return None, ptr


def load_pos_set_single_base(bed_path: str) -> Set[Tuple[str, int]]:
    s: Set[Tuple[str, int]] = set()
    with open(bed_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, _end = line.rstrip("\n").split("\t")[:3]
            chrom = chrom.strip()
            if chrom not in KEEP_CHROMS:
                continue
            s.add((chrom, int(start)))
    return s


def parse_counts_ATCG(counts_field: str) -> Tuple[int, int, int, int]:
    """
    sync counts are A:T:C:G:N:del
    """
    parts = counts_field.split(":")
    while len(parts) < 4:
        parts.append("0")
    A = int(parts[0])
    T = int(parts[1])
    C = int(parts[2])
    G = int(parts[3])
    return A, T, C, G


def ref_alt_mac_depth(ref_base: str, A: int, T: int, C: int, G: int) -> Tuple[int, int]:
    """
    Uses a biallelic approximation:
      - ref count = count(ref)
      - alt count = max(non-ref allele counts)
      - depth_bi = ref + alt
      - mac = min(ref, alt)
    """
    ref_base = ref_base.upper()
    if ref_base not in BASE_TO_IDX:
        raise ValueError("Non-ATCG ref")

    counts = [A, T, C, G]
    ref_i = BASE_TO_IDX[ref_base]
    ref_count = counts[ref_i]

    nonref_counts = [counts[i] for i in range(4) if i != ref_i]
    alt_count = max(nonref_counts) if nonref_counts else 0

    total_nonref = sum(nonref_counts)
    if not ALLOW_MULTI_ALLELIC and total_nonref != alt_count:
        raise ValueError("Multi-allelic site")

    depth_bi = ref_count + alt_count
    if depth_bi <= 0:
        raise ValueError("Zero depth")

    mac = min(ref_count, alt_count)
    return mac, depth_bi


def add_projected_folded(sfs_vec: np.ndarray, mac: int, depth_bi: int, n_proj: int) -> None:
    proj = moments.Numerics._cached_projection(n_proj, depth_bi, mac)
    for i, val in enumerate(proj):
        j = i if i <= (n_proj - i) else (n_proj - i)
        sfs_vec[j] += val


def pop_id_from_filename(fn: str) -> str:
    for suf in [".masked.sync.gz", ".sync.gz", ".gz"]:
        if fn.endswith(suf):
            return fn[: -len(suf)]
    return fn


# -----------------------
# PER-POP COMPUTATION
# -----------------------
def compute_per_gene_sfs_for_pop(
    sync_gz: str,
    pop_id: str,
    genes: Dict[str, List[Tuple[int, int, str]]],
    pos4: Set[Tuple[str, int]],
    pos0: Set[Tuple[str, int]],
    min_depth: int = MIN_DEPTH,   # NEW
):
    gene_sfs_4D = defaultdict(lambda: np.zeros(n_proj + 1, dtype=float))
    gene_sfs_0D = defaultdict(lambda: np.zeros(n_proj + 1, dtype=float))
    gene_L_4D = defaultdict(int)
    gene_L_0D = defaultdict(int)

    ptr = {chrom: 0 for chrom in genes.keys()}

    callable_4D = callable_0D = 0
    poly_sites_4D = poly_sites_0D = 0

    # Optional extra stats (useful to sanity-check min_depth impact)
    not_callable_min_depth = 0

    with gzip.open(sync_gz, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, pos1, ref, counts_field = parts[:4]

            if chrom not in genes:
                continue

            if counts_field.startswith("."):
                continue

            try:
                pos0_coord = int(pos1) - 1
            except ValueError:
                continue

            key = (chrom, pos0_coord)
            is4 = key in pos4
            is0 = key in pos0
            if not (is4 or is0):
                continue

            gid, ptr[chrom] = find_gene_for_pos0(genes[chrom], pos0_coord, ptr[chrom])
            if gid is None:
                continue

            #make "callable" explicit via min_depth on A+T+C+G
            try:
                A, T, C, G = parse_counts_ATCG(counts_field)
            except Exception:
                continue

            called_depth = A + T + C + G
            if called_depth < min_depth:
                not_callable_min_depth += 1
                continue

            # Now compute mac/depth used for projection (and possible multi-allelic filtering)
            try:
                mac, depth_bi = ref_alt_mac_depth(ref, A, T, C, G)
            except Exception:
                continue

            # Count callable length per gene / class
            if is4:
                gene_L_4D[gid] += 1
                callable_4D += 1
            else:
                gene_L_0D[gid] += 1
                callable_0D += 1

            # Monomorphic sites contribute later to bin 0
            if mac == 0:
                continue

            # Polymorphic sites add to projected folded SFS
            if is4:
                poly_sites_4D += 1
                add_projected_folded(gene_sfs_4D[gid], mac, depth_bi, n_proj)
            else:
                poly_sites_0D += 1
                add_projected_folded(gene_sfs_0D[gid], mac, depth_bi, n_proj)

    # fill bin0 per gene
    for gid, L in gene_L_4D.items():
        s = gene_sfs_4D[gid]
        s[0] = max(0.0, float(L) - float(s[1:].sum()))

    for gid, L in gene_L_0D.items():
        s = gene_sfs_0D[gid]
        s[0] = max(0.0, float(L) - float(s[1:].sum()))

    stats = {
        "pop_id": pop_id,
        "min_depth_ATCG": int(min_depth),
        "callable_4D_in_genes": callable_4D,
        "callable_0D_in_genes": callable_0D,
        "poly_4D_sites_in_genes": poly_sites_4D,
        "poly_0D_sites_in_genes": poly_sites_0D,
        "genes_with_4D": len(gene_L_4D),
        "genes_with_0D": len(gene_L_0D),
        "sites_filtered_called_depth_lt_min_depth": not_callable_min_depth,
    }

    return (
        dict(gene_sfs_4D),
        dict(gene_sfs_0D),
        dict(gene_L_4D),
        dict(gene_L_0D),
        stats,
    )


# -----------------------
# WRITE OUTPUTS
# -----------------------
def write_pop_outputs(
    pop_out_dir: str,
    pop_id: str,
    gene_sfs_4D,
    gene_sfs_0D,
    gene_L_4D,
    gene_L_0D,
    stats,
):
    os.makedirs(pop_out_dir, exist_ok=True)

    np.save(os.path.join(pop_out_dir, f"gene_sfs_4D_{pop_id}.npy"), gene_sfs_4D, allow_pickle=True)
    np.save(os.path.join(pop_out_dir, f"gene_sfs_0D_{pop_id}.npy"), gene_sfs_0D, allow_pickle=True)

    with open(os.path.join(pop_out_dir, f"gene_L_4D_{pop_id}.tsv"), "w") as out:
        out.write("gene\tL_4D_callable\n")
        for gid, L in gene_L_4D.items():
            out.write(f"{gid}\t{L}\n")

    with open(os.path.join(pop_out_dir, f"gene_L_0D_{pop_id}.tsv"), "w") as out:
        out.write("gene\tL_0D_callable\n")
        for gid, L in gene_L_0D.items():
            out.write(f"{gid}\t{L}\n")

    with open(os.path.join(pop_out_dir, f"stats_{pop_id}.txt"), "w") as out:
        for k, v in stats.items():
            out.write(f"{k}\t{v}\n")


# -----------------------
# DRIVER
# -----------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading genes:", GENE_BED)
    genes = load_genes_bed(GENE_BED)
    print("  total genes:", sum(len(v) for v in genes.values()))

    print("Loading 4D positions:", FOURFOLD_BED)
    pos4 = load_pos_set_single_base(FOURFOLD_BED)
    print("  4D sites:", len(pos4))

    print("Loading 0D positions:", ZEROFOLD_BED)
    pos0 = load_pos_set_single_base(ZEROFOLD_BED)
    print("  0D sites:", len(pos0))

    sync_files = sorted(fn for fn in os.listdir(SYNC_DIR) if fn.endswith(".masked.sync.gz"))
    print("\nFound", len(sync_files), "masked.sync.gz files")

    for fn in sync_files:
        sync_path = os.path.join(SYNC_DIR, fn)
        pop_id = pop_id_from_filename(fn)
        pop_out_dir = os.path.join(OUT_DIR, pop_id)

        print("\n=== POP:", pop_id, "===")
        gene_sfs_4D, gene_sfs_0D, gene_L_4D, gene_L_0D, stats = compute_per_gene_sfs_for_pop(
            sync_path, pop_id, genes, pos4, pos0, min_depth=MIN_DEPTH
        )

        write_pop_outputs(pop_out_dir, pop_id, gene_sfs_4D, gene_sfs_0D, gene_L_4D, gene_L_0D, stats)

        print("Wrote:", pop_out_dir)
        print("  min_depth_ATCG:", stats["min_depth_ATCG"])
        print("  callable_4D_in_genes:", stats["callable_4D_in_genes"])
        print("  callable_0D_in_genes:", stats["callable_0D_in_genes"])
        print("  poly_4D_sites_in_genes:", stats["poly_4D_sites_in_genes"])
        print("  poly_0D_sites_in_genes:", stats["poly_0D_sites_in_genes"])
        print("  sites_filtered_called_depth_lt_min_depth:", stats["sites_filtered_called_depth_lt_min_depth"])

    print("\nDone. Outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
