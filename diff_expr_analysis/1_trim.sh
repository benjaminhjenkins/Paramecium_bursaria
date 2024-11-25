a=/scratch/ctsII_sequencing/raw_data/220701_A00711_0575_AHKW2YDSX3
c=polyA/

mkdir -p /scratch/ctsII_sequencing/trimmed_reads/$c
cd $a
parallel trim_galore --illumina --paired --fastqc -o /scratch/ctsII_sequencing/trimmed_reads/$c {} {=s/_1/_2/=} ::: *_1.fastq.gz
cd ../