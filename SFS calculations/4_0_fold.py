import gzip

DEGEN_BED = "degenotate_ncbi_out/degeneracy-all-sites.bed"
OUT_4D = "fourfold_sites.bed"
OUT_0D = "zerofold_sites.bed"

with open(DEGEN_BED, "rt") as f, open(OUT_4D, "w") as o4, open(OUT_0D, "w") as o0:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, site_id, deg, *rest = line.rstrip("\n").split("\t")
        if deg == "4":
            o4.write(f"{chrom}\t{start}\t{end}\n")
        elif deg == "0":
            o0.write(f"{chrom}\t{start}\t{end}\n")

print("Wrote", OUT_4D, "and", OUT_0D)
