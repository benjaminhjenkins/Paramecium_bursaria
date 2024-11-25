cd /scratch/ctsII_sequencing/trimmed_reads/polyA

ls *_1_val_1.fq.gz | cut -f 1-3 -d '_' | parallel -j 20 "bowtie2 --very-sensitive-local -x /scratch/ctsII_sequencing/isoseq_assembly/Paramecium_bursaria_186b.mrna-transcripts.minus_CLP_constructonly.fasta -p 4 -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz -S {}_mapped_and_unmapped.sam"

for i in `ls *sam | cut -f 1-3 -d '_'`; do samtools view -bS $i'_mapped_and_unmapped.sam' > $i'_mapped_and_unmapped.bam'; done

ls *mapped_*bam | cut -f 1-3 -d '_' | parallel -j 20 "samtools view -b -f 3 -F 12 {}_mapped_and_unmapped.bam > {}_bothReadsMapped.bam"

ls *both*bam | parallel -j 20 "samtools sort -n -m 100G -@ 4 {} -o {}_sorted.bam"

ls *sorted.bam | cut -f 1-3 -d '_' | parallel -j 20 "samtools fastq -@ 4 {}_bothReadsMapped.bam_sorted.bam -1 {}_host_cleaned_R1.fastq.gz -2 {}_host_cleaned_R2.fastq.gz -0 /dev/null -s /dev/null -n"