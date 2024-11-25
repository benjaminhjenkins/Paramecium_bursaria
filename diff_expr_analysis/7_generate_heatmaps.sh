cd /scratch/ctsII_sequencing/DESeq2_fourway_CLP_cyclo_v1
analyze_diff_expr.pl \
      --matrix /scratch/ctsII_sequencing/salmon_v5/fourway_CLP_cyclo/isoseq_minus_CLPc/Trinity.isoform.TMM.EXPR.matrix \
      --samples /scratch/ctsII_sequencing/trinity/polyA_samples_trimmed.tsv \
      -P 1e-3 -C 0.5