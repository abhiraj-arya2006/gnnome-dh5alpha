"""Annotate DH5α (CP017100.1) repeats.

Inputs : CP017100.1.gb (NCBI, sequence verified identical to CP017100.fasta), nucmer self-alignment
         (show-coords -rclTH -I 85 -L 300), dustmasker intervals.
Outputs: repeats.bed         merged repeat loci (>=300 bp, >=85 % identity to another genomic copy)
         repeat_families.tsv families (copies linked by self-alignment, union-find)
         elements.bed        annotated mobile/structural elements (rRNA, IS/transposase, Rhs, phage, CRISPR, tRNA)
         homopolymers.bed    homopolymer runs >= 7 bp
         lowcomplexity.bed   dustmasker intervals
         summary.json
"""
import sys, re, json
from collections import defaultdict, Counter
from Bio import SeqIO

gb, coords, dust, outdir = sys.argv[1:5]
rec = next(SeqIO.parse(gb, "genbank")); seq = str(rec.seq).upper(); G = len(seq)
chrom = "CP017100.1"

# ---------- annotated elements ----------
def elem_class(f):
    p = (" ".join(f.qualifiers.get("product", [])) + " " + " ".join(f.qualifiers.get("gene", []))).lower()
    if f.type == "rRNA": return "rRNA"
    if f.type == "tRNA": return "tRNA"
    if f.type == "repeat_region": return "CRISPR" if "CRISPR" in str(f.qualifiers.get("rpt_family", "")) else "repeat_region"
    if f.type == "mobile_element": return "IS_transposase"
    if f.type != "CDS": return None
    if "transposase" in p or re.search(r"\bis\d+[a-z]?\b", p) or "insertion element" in p: return "IS_transposase"
    if re.search(r"\brhs", p): return "Rhs"
    if "phage" in p or "integrase" in p or "excisionase" in p: return "phage_related"
    return None
elems = []
for f in rec.features:
    c = elem_class(f)
    if c:
        name = (f.qualifiers.get("product") or f.qualifiers.get("gene") or f.qualifiers.get("rpt_family") or [c])[0]
        elems.append((int(f.location.start), int(f.location.end), c, name.replace("\t", " ")))
elems.sort()

# ---------- self-alignment repeats ----------
hits = []
for l in open(coords):
    s1, e1, s2, e2, l1, l2, idy = l.split("\t")[:7]
    s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
    a2, b2 = min(s2, e2), max(s2, e2)
    if s1 == a2 and e1 == b2: continue                       # trivial diagonal
    hits.append((s1 - 1, e1, a2 - 1, b2, float(idy), "inverted" if s2 > e2 else "direct"))
iv = []
for k, (a1, b1, a2, b2, idy, ori) in enumerate(hits):
    iv.append((a1, b1, k)); iv.append((a2, b2, k))
iv.sort()
loci = []                                                    # merge overlapping intervals into loci
for a, b, k in iv:
    if loci and a <= loci[-1][1]:
        loci[-1][1] = max(loci[-1][1], b); loci[-1][2].add(k)
    else: loci.append([a, b, {k}])
hit2loci = defaultdict(set)
for i, (a, b, ks) in enumerate(loci):
    for k in ks: hit2loci[k].add(i)
parent = list(range(len(loci)))
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for k, ls in hit2loci.items():
    ls = list(ls)
    for o in ls[1:]: parent[find(o)] = find(ls[0])
fam_members = defaultdict(list)
for i in range(len(loci)): fam_members[find(i)].append(i)

def overlap_classes(a, b):
    c = Counter()
    for ea, eb, cl, nm in elems:
        if ea < b and eb > a: c[cl] += min(b, eb) - max(a, ea)
    return c
locus_info = []
fam_rows = []
fam_id = {}
for n, (root, mem) in enumerate(sorted(fam_members.items(), key=lambda x: -len(x[1]))):
    cls = Counter(); maxid = 0; oris = set(); lens = []
    for i in mem:
        a, b, ks = loci[i]; cls.update(overlap_classes(a, b)); lens.append(b - a)
        for k in ks: maxid = max(maxid, hits[k][4]); oris.add(hits[k][5])
    # family class priority: rRNA > IS > Rhs > phage > CRISPR > tRNA > unannotated
    order = ["rRNA", "IS_transposase", "Rhs", "phage_related", "CRISPR", "tRNA", "repeat_region"]
    fclass = next((o for o in order if cls.get(o, 0) >= 100), "unannotated")
    fid = f"F{n+1:03d}"
    for i in mem: fam_id[i] = (fid, fclass, len(mem), maxid)
    fam_rows.append((fid, fclass, len(mem), min(lens), max(lens), sum(lens), maxid, ",".join(sorted(oris)),
                     ";".join(f"{loci[i][0]}-{loci[i][1]}" for i in sorted(mem, key=lambda i: loci[i][0]))))
with open(f"{outdir}/repeats.bed", "w") as w:
    for i, (a, b, ks) in enumerate(loci):
        fid, fc, cn, mi = fam_id[i]
        w.write(f"{chrom}\t{a}\t{b}\tL{i+1:03d}\t{fid}\t{fc}\t{cn}\t{mi:.2f}\n")
with open(f"{outdir}/repeat_families.tsv", "w") as w:
    w.write("family\tclass\tcopies\tmin_len\tmax_len\ttotal_bp\tmax_identity\torientations\tloci\n")
    for r in fam_rows: w.write("\t".join(str(x) for x in r) + "\n")
with open(f"{outdir}/elements.bed", "w") as w:
    for a, b, c, nm in elems: w.write(f"{chrom}\t{a}\t{b}\t{c}\t{nm}\n")
# ---------- homopolymers and low complexity ----------
hp = [(m.start(), m.end(), m.group()[0]) for m in re.finditer(r"A{7,}|C{7,}|G{7,}|T{7,}", seq)]
with open(f"{outdir}/homopolymers.bed", "w") as w:
    for a, b, base in hp: w.write(f"{chrom}\t{a}\t{b}\t{base}{b-a}\n")
dl = []
for l in open(dust):
    if l.startswith(">"): continue
    a, b = l.strip().split(" - "); dl.append((int(a), int(b) + 1))
with open(f"{outdir}/lowcomplexity.bed", "w") as w:
    for a, b in dl: w.write(f"{chrom}\t{a}\t{b}\tdust\n")
rep_bp = sum(b - a for a, b, _ in loci)
S = dict(genome_len=G, self_alignment_hits=len(hits), repeat_loci=len(loci), repeat_families=len(fam_rows),
         repeat_bp=rep_bp, repeat_frac=rep_bp / G,
         family_classes=dict(Counter(r[1] for r in fam_rows)),
         repeat_bp_by_class={c: sum(loci[i][1] - loci[i][0] for i in range(len(loci)) if fam_id[i][1] == c) for c in set(r[1] for r in fam_rows)},
         elements_by_class=dict(Counter(c for _, _, c, _ in elems)),
         homopolymers_ge7=len(hp), homopolymers_ge8=sum(b - a >= 8 for a, b, _ in hp), lowcomplexity_intervals=len(dl),
         lowcomplexity_bp=sum(b - a for a, b in dl))
json.dump(S, open(f"{outdir}/summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
for r in fam_rows[:30]: print("\t".join(str(x) for x in r[:8]))
