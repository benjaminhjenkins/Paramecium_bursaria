#Replace "seq" with sequence ID

awk '/^>/{f=!d[$1];d[$1]=1}f' seq.fasta > seq_d.fasta
seqtk seq -L 100 seq_d.fasta > seq_d_l100.fasta

mafft seq_d_l100.fasta > seq_d_l100_aligned.fasta
trimal -in seq_d_l100_aligned.fasta -out seq_d_l100_masked_0.9.fasta -gt 0.9
iqtree -s seq_d_l100_masked_0.9.fasta -nt AUTO -bb 1000 -alrt 1000
