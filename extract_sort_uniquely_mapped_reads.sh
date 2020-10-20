BAM_FILES=$(ls *.bam |cut -d "." -f1)


for i in $BAM_FILES; do samtools view -H ${i}.bam > ${i}_header.sam; samtools view ${i}.bam |grep -w "NH:i:1" | cat ${i}_header.sam - |samtools view -@ 10 -Sb - > ${i}_unique.bam; samtools sort -@ 10 ${i}_unique.bam > ${i}_unique_sorted.bam; rm ${i}.sam; rm ${i}_header.sam; rm ${i}.bam; rm ${i}_unique.bam; done

