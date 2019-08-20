#!/bin/bash
# reference_base_assembly_pipeline
sample_prefix=SP0001
read1='$sample_prefix'_R1_val_1.fq.gz'
read2='$sample_prefix'_R2_val_2.fq.gz'
fasta_file=1045684451.fasta
#bwa mapping
bwa mem $fasta_file $read1 $read2 > $sample_prefix'.sam'
#samtools sort
samtools sort -O bam -T temp1 $sample_prefix'.sam' >| $sample_prefix'.bam'
#samtools index
samtools index $sample_prefix'.bam'
#samtools mpileup
samtools mpileup -f 1045684451.fasta -gu $sample_prefix'.bam' | bcftools call -c -O b -o
$sample_prefix'.raw.bcf'
#convert file to fastq format
bcftools view -O v $sample_prefix'.raw.bcf' | vcfutils.pl vcf2fq > $sample_prefix'.fastq'
#convert fastq to fasta
python3 convert_fastq_to_fasta.py -q $sample_prefix'.fastq' -a $sample_prefix'.fasta'
