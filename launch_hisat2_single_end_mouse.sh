HISAT2_INDEXES=/home/ESPRESSO/antonio/rna-seq_tools/indexes/grcm38
FASTQ_FILES=$(ls *.fastq.gz |cut -d "." -f1)

for i in $FASTQ_FILES; do hisat2 -q -p 20 -x $HISAT2_INDEXES/genome -U ${i}.fastq.gz | samtools view -Sbh > ${i}.bam; done
