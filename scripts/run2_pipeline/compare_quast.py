import sys, csv
keys = ["# contigs","Total length","Reference length","Genome fraction (%)","Duplication ratio","# misassemblies","# local misassemblies","Unaligned length","NGA50","# mismatches per 100 kbp","# indels per 100 kbp","Largest alignment"]
cols = []
for arg in sys.argv[1:]:
    label, path = arg.split("=", 1)
    d = {r[0]: r[1] for r in csv.reader(open(path), delimiter="\t") if len(r) > 1}
    cols.append((label, d))
w = max(len(k) for k in keys) + 2
print("Metric".ljust(w) + "".join(l.rjust(22) for l, _ in cols))
for k in keys:
    print(k.ljust(w) + "".join(d.get(k, "-").rjust(22) for _, d in cols))
