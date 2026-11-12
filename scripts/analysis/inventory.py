import hashlib, os, sys
B = "/home/abhiraj_arya/research/"
rows = [
 # label, run, status, assembly, quast_dir, notes
 ("R1_Raven_p0","Run1","VALID","GNNome/data/CP017100_assembly.fasta","GNNome/data/classical_raven_quast_results","Raven fork -p0, 60x HiFi-like reads"),
 ("R1_GNNome","Run1","VALID","GNNome/data/raven/assembly/0_assembly.fasta","GNNome/data/raven/quast_results","GNNome on Run1 graph"),
 ("R2old_Raven_p2_conda","Run2-old","VALID (reproduced byte-identical)","GNNome/data/dh5alpha/raven/output/dh5alpha_asm.fasta","GNNome/data/dh5alpha/classical_raven_quast_results","Conda Raven -p2; == R2_Raven_internal_p2"),
 ("R2old_GNNome_K12graph","Run2-old","INVALID (K-12 graph)","GNNome/data/dh5alpha/raven/assembly/0_assembly.fasta","GNNome/data/dh5alpha/quast_gnnome_raven","consumed stale K-12 graph_1.gfa"),
 ("R2_Raven_p0","Run2","VALID","run2_corrected/02_raven_p0/raven_p0_assembly.fasta","run2_corrected/05_quast/minid80_raven_p0","default QUAST cannot align (<95% id); minid80 supplementary"),
 ("R2_GNNome","Run2","VALID","run2_corrected/04_gnnome/raven/assembly/0_assembly.fasta","run2_corrected/05_quast/minid80_gnnome_p0","as above"),
 ("R2_Raven_Racon2","Run2","VALID","run2_corrected/07_polish/raven_racon2/round2.fasta","run2_corrected/07_polish/quast/raven_racon2","external minimap2+Racon x2"),
 ("R2_GNNome_Racon2","Run2","VALID","run2_corrected/07_polish/gnnome_racon2/round2.fasta","run2_corrected/07_polish/quast/gnnome_racon2","external minimap2+Racon x2"),
 ("R2_Raven_internal_p2","Run2","VALID","run2_corrected/07_polish/raven_internal_p2/raven_p2_assembly.fasta","run2_corrected/07_polish/quast/raven_internal_p2","Raven fork -p2"),
]
def fasta_len(p):
    n = c = 0
    for l in open(p):
        if l.startswith(">"): c += 1
        else: n += len(l.strip())
    return c, n
print("label\trun\tstatus\tcontigs\tlength\tsha256\tassembly\tquast_report\tquast_present\tnotes")
for lab, run, st, a, q, note in rows:
    c, n = fasta_len(B + a)
    h = hashlib.sha256(open(B + a, "rb").read()).hexdigest()
    print(f"{lab}\t{run}\t{st}\t{c}\t{n}\t{h}\t{a}\t{q}\t{os.path.isfile(B+q+'/report.tsv')}\t{note}")
