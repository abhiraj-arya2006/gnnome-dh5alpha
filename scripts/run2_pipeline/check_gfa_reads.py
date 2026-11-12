import sys, gzip, random
gfa, *readfiles = sys.argv[1:]
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
segs = {}
for line in open(gfa):
    if line.startswith("S\t"):
        p = line.split("\t", 3); segs[p[1]] = p[2]
random.seed(42)
sample = dict(random.sample(sorted(segs.items()), 200))
def seqs(path):
    f = (gzip.open if path.endswith(".gz") else open)(path, "rt")
    if f.read(1) == "@":
        f.seek(0)
        for h in f:
            s = f.readline().strip(); f.readline(); f.readline()
            yield h[1:].split()[0], s
    else:
        f.seek(0); name, buf = None, []
        for line in f:
            if line.startswith(">"):
                if name: yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else: buf.append(line.strip())
        if name: yield name, "".join(buf)
print(f"{gfa}: {len(segs)} segments, testing 200 random (seed 42)")
for rf in readfiles:
    found = {n: s for n, s in seqs(rf) if n in sample}
    hit = sum(1 for n, v in sample.items() if n in found and (v in found[n] or v in rc(found[n])))
    print(f"  vs {rf}: same read name {len(found)}/200, identical sequence {hit}/200")
