"""Recalculate structural breakpoints of each assembly against CP017100.1.

Two independent sources per assembly:
  quast    : consecutive aligned blocks in QUAST's all_alignments_*.tsv and the label QUAST put between them
  minimap2 : minimap2 -cx asm20 (primary + supplementary), blocks sorted by contig coordinate; each junction is
             classified from the reference jump vs. contig gap:
               inversion   strand changes
               relocation  |ref_jump - contig_gap| > 1000 bp or the reference goes backwards by > 1000 bp
               local       200 < |ref_jump - contig_gap| <= 1000 bp
               origin      junction joins the reference end to the reference start (circular genome, linear simulation)
             Junctions with a discrepancy <= 200 bp are alignment fragmentation, not breakpoints (QUAST's definition).
usage: breakpoints.py <label> <assembly.fasta> <quast_dir or -> <out.tsv>   (needs REF, MM env vars)
"""
import sys, os, subprocess, glob
label, asm, qdir, out = sys.argv[1:5]
REF, MM = os.environ["REF"], os.environ["MM"]
G = 4583637
rows = []

def classify(b1, b2):
    # b = (qs, qe, strand, rs, re)
    if b1[2] != b2[2]: return "inversion", None, None
    qgap = b2[0] - b1[1]
    if b1[2] == "+": rjump = b2[3] - b1[4]
    else: rjump = b1[3] - b2[4]
    d = rjump - qgap
    if (b1[2] == "+" and b1[4] > G - 5000 and b2[3] < 5000) or (b1[2] == "-" and b1[3] < 5000 and b2[4] > G - 5000):
        return "origin", rjump, qgap
    if abs(d) > 1000 or rjump < -1000: return "relocation", rjump, qgap
    if abs(d) > 200: return "local", rjump, qgap
    return None, rjump, qgap

def emit(src, blocks):
    blocks = sorted(blocks, key=lambda b: (b[5], b[0]))
    for b1, b2 in zip(blocks, blocks[1:]):
        if b1[5] != b2[5]: continue
        t, rj, qg = classify(b1, b2)
        if t:
            rows.append((label, src, t, b1[4] if b1[2] == "+" else b1[3], b2[3] if b2[2] == "+" else b2[4], b1[1], b2[0],
                         rj, qg, None if rj is None else rj - qg, f"{b1[6]:.4f}/{b2[6]:.4f}"))
    return len(blocks)
def mm_blocks(extra):
    paf = subprocess.run([MM, "-t", "2", "-cx", "asm20", "--secondary=no"] + extra + [REF, asm], capture_output=True, text=True).stdout
    out = []
    for l in paf.splitlines():
        f = l.split("\t")
        if int(f[11]) < 1 or int(f[10]) < 1000: continue
        out.append((int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8]), f[0], int(f[9]) / int(f[10])))
    return out
nb = {}
nb["minimap2"] = emit("minimap2", mm_blocks([]))
nb["minimap2_nolj"] = emit("minimap2_nolj", mm_blocks(["--no-long-join", "-r", "200"]))
if os.environ.get("NUCMER") == "1":
    tmp = os.path.join(os.path.dirname(out) or ".", "nucmer_" + label)
    subprocess.run(["nucmer", "-p", tmp, REF, asm], capture_output=True)
    subprocess.run(f"delta-filter -1 {tmp}.delta > {tmp}.1delta", shell=True)
    co = subprocess.run(["show-coords", "-THrcl", tmp + ".1delta"], capture_output=True, text=True).stdout
    bl = []
    for l in co.splitlines():
        f = l.split("\t"); s1, e1, s2, e2, l1, l2, idy = int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5]), float(f[6])
        if l2 < 1000: continue
        st = "+" if s2 < e2 else "-"
        bl.append((min(s2, e2) - 1, max(s2, e2), st, s1 - 1, e1, f[-1], idy / 100))
    nb["nucmer"] = emit("nucmer", bl)
    for ext in (".delta", ".1delta"):
        if os.path.exists(tmp + ext): os.remove(tmp + ext)
# ---- QUAST ----
if qdir != "-":
    f = glob.glob(f"{qdir}/contigs_reports/all_alignments_*.tsv")
    if f:
        lines = [l.rstrip("\n") for l in open(f[0])]
        prev = None; between = None
        for l in lines:
            p = l.split("\t")
            if p[0].isdigit():
                s1, e1, s2, e2 = map(int, p[:4]); cur = (s1, e1, s2, e2)
                if prev and between:
                    rows.append((label, "quast", between, prev[1], cur[0], prev[3], cur[2], cur[0] - prev[1],
                                 cur[2] - prev[3], (cur[0] - prev[1]) - (cur[2] - prev[3]), p[6]))
                prev, between = cur, None
            elif p[0] in ("CONTIG", "S1"): prev, between = (prev if p[0] == "S1" else None), None
            else: between = l.strip().replace(" ", "_")
with open(out, "w") as w:
    w.write("assembly\tsource\ttype\tref_left\tref_right\tcontig_left\tcontig_right\tref_jump\tcontig_gap\tdiscrepancy\tidentity\n")
    for r in rows: w.write("\t".join("" if x is None else str(x) for x in r) + "\n")
from collections import Counter
c = Counter((r[1], r[2]) for r in rows)
print(label, dict(c), "blocks>=1kb:", nb)
