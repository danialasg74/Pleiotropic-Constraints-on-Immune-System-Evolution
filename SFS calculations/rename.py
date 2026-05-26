import re

FASTA = "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
IN_BEDS = ["zerofold_sites.bed", "fourfold_sites.bed"]
KEEP = {"2L","2R","3L","3R","X","4","Y"}  # DEST names
OUT_SUFFIX = ".destnames.bed"

# Build accession -> DEST chrom mapping from FASTA headers
acc_to_dest = {}

with open(FASTA, "rt") as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        acc = header.split()[0]  # e.g. NC_004354.4

        # match "chromosome X" or "chromosome 2L" etc
        m = re.search(r"\bchromosome\s+([0-9A-Za-z]+)\b", header)
        if not m:
            continue
        chrom = m.group(1)

        # normalize names to DEST style
        if chrom in KEEP:
            acc_to_dest[acc] = chrom
        elif chrom in {"M", "MT", "mitochondrion"}:
            
            continue

print("Example mappings:", list(acc_to_dest.items())[:5])
print("Total mapped accessions:", len(acc_to_dest))

for bed in IN_BEDS:
    out_bed = bed.replace(".bed", OUT_SUFFIX)
    kept = 0
    skipped = 0

    with open(bed, "rt") as inp, open(out_bed, "wt") as out:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            new = acc_to_dest.get(chrom)
            if new is None:
                skipped += 1
                continue
            out.write(f"{new}\t{start}\t{end}\n")
            kept += 1

    print(f"{bed} -> {out_bed} | kept {kept:,} | skipped {skipped:,}")
