cd /scratch/ctsII_sequencing/salmon_v5/fourway_CLP_cyclo/isoseq_minus_CLPc
align_and_estimate_abundance.pl --seqType fq --samples_file /scratch/ctsII_sequencing/trinity/polyA_samples_trimmed.tsv \
--transcripts /scratch/ctsII_sequencing/isoseq_assembly/Paramecium_bursaria_186b.mrna-transcripts.minus_CLP_constructonly.fasta \
--est_method salmon --trinity_mode --prep_reference