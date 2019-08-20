#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) give reads1 name, like "/path of reads directory/*R1.fastq.gz";
(2) give path of blastn database; here it is "~/softwares/ncbi_database/virus_genomes/all_virus_genome2.fna";
(3) give path of blastx database; here it is "~/softwares/ncbi_database/virus_proteome/viral_refseq.fa.dmnd";
(4) give k-mer size like 21,31,41,51,61,71,81,91,101,121 (if read length >=250) or 21,31,41,51,61,71 (if read length upto 150bp);
')

# Check for the right number of arguments. Assign them to variables.

NumArgsExpected=4
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi

echo
echo =========================================================================
echo Script for Metagenomic analysis 
echo Author: Vikash K. Singh
echo Affliation: University of Minnesota - College of veterinary and population medicine
echo email: vsingh@umn.edu
echo =========================================================================
echo
sleep 5s

WorkDir=`pwd`
reads1=$1
blastn_db=$2
blastx_db=$3

for i in $reads1
do
    filterbytrimmomatic.pl $i ${i/R1/R2}
done

find . -type f -name "*fastq" | parallel pigz -p 50 {}

for i in *R1_paired_trimmed.fastq.gz 
do 
    metaspades.py \
    --only-assembler \
    -t 50 \
    -m 200 \
    -k $4 \
    -o `basename ${i/_R1*/.ASM_out}` \
    --pe1-1 $i \
    --pe1-2 ${i/R1/R2} \
    --pe1-s ${i/paired/single} \
    --pe1-s ${i/R1_paired/R2_single} &
done

wait

for i in $WorkDir/*ASM_out
do
    perl -e ' $count=0; $len=0; while(<>) { s/\r?\n//; s/\t/ /g; if (s/^>//) { if ($. != 1) { print "\n" } s/ |$/\t/; $count++; $_ .= "\t"; } else { s/ //g; $len += length($_) } print $_; } print "\n"; warn "\nConverted $count FASTA records in $. lines to tabular format\nTotal sequence length: $len\n\n"; ' $i/scaffolds.fasta > `basename ${i/ASM*/fasta.tab}`
done 

for i in $WorkDir/*ASM_out
do 
    blastn -query $i/scaffolds.fasta \
    -db $2 \
    -out `basename ${i/out/blastn}` \
    -outfmt "6 qseqid qlen qcovs sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 1e-06 \
    -max_target_seqs 1 \
    -num_threads 30
done

for i in $WorkDir/*ASM_out
do 
    diamond blastx \
    -d $3 \
    -q $i/scaffolds.fasta \
    -f 6 qseqid qlen qcovhsp stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
    --sensitive \
    --no-auto-append \
    --top 10 \
    --out ${i/out/blastx}
done

for i in $WorkDir/*blastn; do cat ~/softwares/myscripts/blastn_header $i > $i.2;done
for i in $WorkDir/*blastx; do cat ~/softwares/myscripts/blastx_header $i > $i.2;done
rm *blastn
rm *blastx
for i in $WorkDir/*blastn.2;do mv $i ${i/blastn*/blastn};done
for i in $WorkDir/*blastx.2;do mv $i ${i/blastx*/blastx};done

ssconvert --merge-to="${PWD##*/}_blastn.xls" *blastn
ssconvert --merge-to="${PWD##*/}_blastx.xls" *blastx
ssconvert --merge-to="${PWD##*/}_fasta.xls" *tab
