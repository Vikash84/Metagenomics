#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) give reads1 name, like "/path of reads directory/*R1.fastq.gz";
(2) give path of ref fasta file; like "/path of ref fasta/ref.fasta";
(3) --local or --sensitive;
')

# Check for the right number of arguments. Assign them to variables.

NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi

echo
echo =========================================================================
echo Script for reference-based assembly 
echo Author: Vikash K. Singh
echo Affliation: University of Minnesota - College of veterinary and population medicine
echo email: vsingh@umn.edu
echo =========================================================================
echo
sleep 5s


WorkDir=`pwd`
reads1=$1
ref_fasta=$2

for i in $reads1
do
    filterbytrimmomatic.pl $i ${i/R1/R2}
done



conda activate /home/vsingh/miniconda3/envs/assembly



for i in $ref_fasta
do 
    samtools faidx $i
done

for i in $ref_fasta
do
    bowtie2-build $i ${i/.fasta/}
done

for i in *R1_paired_trimmed.fastq
do
    bowtie2 \
    $3 \
    -p 50 \
    -x ${ref_fasta/.fasta} \
    -1 $i \
    -2 ${i/R1/R2} \
    -S `basename ${i/_R1*/.sam}`
done

for i in *.sam
do 
    samtools view \
    -bS $i | samtools sort - -o ${i/sam/bam}
done

for i in *.bam;do
    samtools index $i
done

for i in *.bam
do 
    [ ! -f ${i/bam/consensus.fa} ] \
    && samtools mpileup --max-depth 10000 \
    -ud 1000 \
    -f $ref_fasta $i | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - > ${i/bam/consensus.fa}
done

# bcftools mpileup --max-depth 10000 -Ou -ABf vaccine.fasta D19-028511-1_VRIS_S1.bam | bcftools call -c -p 0.05 -P 1.1e-5 | vcfutils.pl vcf2fq | seqtk seq -a - > mapped.fasta

for i in *bam;do /home/vsingh/miniconda3/envs/assembly/bin/java -Xmx100G -jar /home/vsingh/miniconda3/envs/assembly/share/pilon-1.23-0/pilon-1.23.jar --genome $2 --fix all --changes --frags $i --threads 50 --output ${i/bam/Pilon_round2};done
