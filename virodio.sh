#!/bin/bash
# Made by Thomas Arn Hansen 2016/2017
# tool for explorative analysis of metagenomics sequence data, mainly for virus discovery
# Virus discovery using Illumina sequencing can be dependent on few reads or small contigs with low homology to known references. This is a wrapper that uses different availble search tools (Kaiju, kraken, blastn and diamond) for explorative virus discovery. The translated searches are usually used in virus discovery to give higher sensitivity.

#variables of tools
trim_galore=/path/to/trim_galore/
bowtie2=/path/to/bowtie2-2.2.9/
kraken=/path/to/kraken/
krona=/path/to/Krona/KronaTools/scripts/
kaiju=/path/to/kaiju/bin/
megahit=/path/to/megahit/
blastn=/path/to/ncbi-blast-2.3.0+/bin/
$diamond=/path/to/diamond

#database Variables:
# bowtie2 database of host organism
depletion_db=/path/to/db/Bowtie2/hg/hg19/hg19
kraken_db=/path/to/kraken/KrakenDB/
kaiju_nodes=/path/to/kaiju/kaijudb/nodes.dmp
kaiju_db=/path/to/kaiju/kaijudb/kaiju_db.fmi
blastn_db=/path/to/ncbi-blast-2.3.0+/db/nt
# fasta for diamond database get it here : ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond_db=/path/to/Diamond/nr.nowhite.dmnd

###################################
## Read in arguments
###############################
if [ $# -ne 3 ]; 
then
echo -e "Wrong number of input arguments\nUsage: sh virodio.sh pair1=/dir/to/inputfile_pair1.txt pair2=/dir/to/inputfile_pair2.txt outdir=/dir/to/outdir";
exit
fi

for var in "$@";
do

eval  $var
done


if [ ! -f $pair1  -o ! -f $pair2 -o ! -d $outdir ] || [ $pair1 = $pair2 ]; then 
echo -e "File not found or pair1 is identical to pair2! or files do not exist  \nUsage: sh virodio.sh outdir=/dir/to/output pair1=/dir/to/pair1.fq pair2=/dir/to/pair2.fq index=/dir/to/index/file genomes=/home/arn/data/genomes/hg19/hg19"
exit
fi

###################################
## Run trim galore -- QC and Fastqc
###############################
mkdir $outdir/0_trim_galore/
$trim_galore/trim_galore --paired --no_report_file --dont_gzip -q 30 --retain_unpaired $pair1 $pair2 -o $outdir/0_trim_galore/
fastqc $pair1 $pair2 -o $outdir/0_trim_galore/

###################################
########## Host depletion
###############################
mkdir $outdir/1_human_depleted/
mkdir $outdir/1_human_depleted/human/
mkdir $outdir/1_human_depleted/NotHuman/

bsname=`basename $outdir/0_trim_galore/*val_1*\.fq _R1_001_val_1.fq`

pth_human=`echo $outdir/1_human_depleted/human/`
pth_Nothuman=`echo $outdir/1_human_depleted/NotHuman/`

fil1=`echo $outdir/0_trim_galore/$bsname\_R1_001_val_1.fq`
fil2=`echo $outdir/0_trim_galore/$bsname\_R2_001_val_2.fq`
fil3=`echo $outdir/0_trim_galore/$bsname\_R1_001_unpaired_1.fq`
fil4=`echo $outdir/0_trim_galore/$bsname\_R2_001_unpaired_2.fq`

$bowtie2/bowtie2 --no-unal --end-to-end -p 4 -q -k 1 -x $depletion_db -1 $fil1 -2 $fil2 -U $fil3,$fil4 --un-conc $pth_Nothuman/$bsname.fastq --un $pth_Nothuman/$bsname.singelton.fastq -S $pth_human/$bsname.sam


###################################
## Run Kraken for classification of reads
###############################
mkdir $outdir/2_kraken/
file1=`ls $pth_Nothuman/$bsname.1.fastq`
file2=`ls $pth_Nothuman/$bsname.2.fastq`

nice -n 19 $kraken/./kraken --db $kraken_db --threads 10 --fastq-input --only-classified-output --output $outdir/2_kraken/output.txt --paired  $file1 $file2

cut -f2,3 $outdir/2_kraken/output.txt > $outdir/2_kraken/krona.int 
perl $krona/ImportTaxonomy.pl $outdir/2_kraken/krona.int  -o $outdir/2_kraken/krona.out.html


###################################
## Run kaiju for classification of reads
###############################
mkdir $outdir/3_kaiju/
$kaiju/./kaiju -t $kaiju_nodes -f $kaiju_db -i $file1 -j $file2 -o $outdir/3_kaiju/$bsname -z 10 -a greedy -e 5

cut -f2,3 $outdir/3_kaiju/$bsname > $outdir/3_kaiju/$bsname.krona.int 
perl $krona/ImportTaxonomy.pl $outdir/3_kaiju/$bsname.krona.int -o $outdir/3_kaiju/$bsname.krona.out.html 

##############################################################
############################## Do megahit! ###################
##############################################################

bsname=`basename $outdir/0_trim_galore/*val_1*\.fq _R1_001_val_1.fq`
echo $bsname

pth_human=`echo $outdir/1_human_depleted/human/`
pth_Nothuman=`echo $outdir/1_human_depleted/NotHuman/`

mkdir $outdir/4_megahit/

fil1=`echo $pth_Nothuman/$bsname.1.fastq`
fil2=`echo $pth_Nothuman/$bsname.2.fastq`
fil3=`echo $pth_Nothuman/$bsname.singelton.fastq`

echo $fil1
echo $fil2
echo $fil3

$megahit./megahit --presets meta-sensitive --min-contig-len 100 -1 $fil1 -2 $fil2 -r $fil3 -o $outdir/4_megahit/$bsname


###################################
## Do BLASTn fo classification
###############################
mkdir $outdir/5_blast/
$blastn./blastn -outfmt 6 -num_threads 16 -query $outdir/4_megahit/$bsname/final.contigs.fa -out $outdir/5_blast/$bsname.tsv -db $blastn_db

perl $krona/ImportBLAST.pl $outdir/5_blast/$bsname.tsv -o $outdir/5_blast/$bsname.test.krona.html


###################################
## Do Diamond fo classification
###############################

mkdir $outdir/6_Diamond/
mkdir $outdir/tmp/

nice -n 19 $diamond/diamond blastx --sensitive -d $diamond_db -q $outdir/4_megahit/$bsname/final.contigs.fa -a $outdir/6_Diamond/$bsname.contigs -t $outdir/tmp/
nice -n 19 $diamond/diamond view -f tab -a $outdir/6_Diamond/$bsname.contigs.daa  -o $outdir/6_Diamond/$bsname.contigs.m8

perl $krona/ImportBLAST.pl $outdir/6_Diamond/$bsname.contigs.m8 -o $outdir/6_Diamond/$bsname.test.krona.html



