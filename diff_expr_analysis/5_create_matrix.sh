cd /scratch/ctsII_sequencing/salmon_v5/fourway_CLP_cyclo/isoseq_minus_CLPc
find S* C* -name "quant.sf" | tee quant_files.list
abundance_estimates_to_matrix.pl --est_method salmon \
--out_prefix Trinity --name_sample_by_basedir \
--quant_files quant_files.list \
--gene_trans_map none