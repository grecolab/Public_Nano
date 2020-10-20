HISAT2_INDEXES=/home/ESPRESSO/antonio/rna-seq_tools/indexes/grcm38
FASTQ_FILES_R1=$(ls *_1.fastq.gz |cut -d "_" -f1)
FASTQ_FILES_R2=$(ls *_2.fastq.gz |cut -d "_" -f1)

for i in $FASTQ_FILES_R1; do hisat2 -q -p 20 -x $HISAT2_INDEXES/genome -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz | samtools view -Sbh > ${i}.bam; done
