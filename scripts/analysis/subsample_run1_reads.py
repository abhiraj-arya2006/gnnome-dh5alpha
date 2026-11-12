"""Seeded random subsample of the Run 1 reads with the matching MAF blocks (exact truth kept).
usage: subsample_run1_reads.py <fraction> <seed> <out_prefix>"""
import sys, gzip, random
frac, seed, out = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
FQ = "/home/abhiraj_arya/research/GNNome/data/CP017100_sim.fastq"; MAF = "/home/abhiraj_arya/research/GNNome/data/CP017100_sim_0001.maf.gz"
ids = [l[1:].split()[0] for i, l in enumerate(open(FQ)) if i % 4 == 0]
rng = random.Random(seed); keep = set(r for r in ids if rng.random() < frac)
with open(FQ) as f, gzip.open(out + "_0001.fq.gz", "wt") as w:
    for h in f:
        s, p, q = f.readline(), f.readline(), f.readline()
        if h[1:].split()[0] in keep: w.write(h + s + p + q)
with gzip.open(MAF, "rt") as f, gzip.open(out + "_0001.maf.gz", "wt") as w:
    block = []
    def flush(b):
        s_lines = [x for x in b if x.startswith("s ")]
        if len(s_lines) == 2 and s_lines[1].split()[1] in keep: w.write("".join(b))
    for l in f:
        if l.startswith("a"):
            if block: flush(block)
            block = [l]
        elif block: block.append(l)
    if block: flush(block)
print(f"kept {len(keep)} of {len(ids)} reads")
