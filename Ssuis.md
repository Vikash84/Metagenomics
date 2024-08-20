# H.parasuis Analysis

### sample sheet preparation
```bash
for i in *R1_001.fastq.gz;do echo -e "${i/_S*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples.tab;done

for i in *fasta;do echo -e "${i/.fasta*}\tGENOME\t`pwd`/$i" >> samples.tab;done
```
### Downloading refence sequenece from NCBI

### Running pipeline
```bash
conda deactivate
PERL5LIB="";
conda activate nullarbor-new

rm ~/miniconda3/envs/nullarbor-new/db/vfdb/sequences*
cp ~/miniconda3/envs/nullarbor-new/db/vfdb/vfdb_all/sequences* ~/miniconda3/envs/nullarbor-new/db/vfdb

PERL5LIB="";nullarbor.pl --cpus 50 --run --trim --mlst hparasuis --treebuilder iqtree_slow --taxoner centrifuge --name Hparasuis_analysis --ref ref.fasta --input samples.tab --outdir results
```

### Serotyping
```bash
blastn -task blastn-short -query Hparasuis_sero_primers.fasta -db HPS.all_merged.fasta -out HPS.all_merged.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50

sort -k3,3nr -k15,15nr HPS.all_merged.blastn | awk -F"\t" '!seen[$1]++' | less


alternate method
conda deactivate
PERL5LIB="";
conda activate GparasuisSero-env
for i in *fasta;do HpsuisSero.sh -i $i -x fasta -o . -s ${i/.fasta};done
```
### vtaA profile
```bash
blastn -task blastn-short -query vtaA_profile.fasta -db scaffolds.fasta -out vtaA_profile.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50

sort -k3,3nr -k15,15nr vtaA_profile.blastn | awk -F"\t" '!seen[$1]++' | less

Or with complete vtaA sequences
blastn -task blastn-short -query vtaA-complete.fasta -db D22-004701-2-repeat.fasta -out HPS.all_merged.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50
```

```bash
for i in *R1*fastq;do tanoti -P 50 -p 1 -r vtaA-complete.fasta -i $i ${i/R1/R2} -o ${i/_S*/.sam};done
for i in *sam;do SAM2CONSENSUS -i $i -o ${i/sam/fasta};done
```



# SSuis Analysis
```bash
for i in *R1_001.fastq.gz;do echo -e "${i/_S*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples_"$( date +"%Y-%m-%d" )".tab;done

nullarbor.pl --cpus 50 --run --mlst ssuis --treebuilder iqtree_slow --taxoner centrifuge --name Ssuis_analysis-"$( date +"%Y-%m-%d" )" --ref Ssuis_ref.fa --input samples_"$( date +"%Y-%m-%d" )".tab --outdir Ssuis_analysis-"$( date +"%Y-%m-%d" )" --assembler spades
```
### Ssuis serotyping
```bash
conda deactivate
PERL5LIB="";
conda activate ssuissero-env
for i in *.fasta;do SsuisSero.sh -s ${i/.fasta} -i $i -o . -x fasta;done
cat *tsv | grep -v Sample | less
```
or

```bash
git clone https://github.com/streplab/SsuisSerotyping_pipeline.git
conda activate srst2-env
for i in *L001_R1_001.fastq;do mv $i ${i/L001_R1_001.fastq/R1_001.fastq};done
for i in *L001_R2_001.fastq;do mv $i ${i/L001_R2_001.fastq/R2_001.fastq};done
```
#
```bash
./Ssuis_serotypingPipeline.pl --fastq_directory `pwd`/data --forward _R1_001 --reverse _R2_001 --ends pe
```

# tormes-1.3.0 E.coli and salmonella analysis
```bash
conda deactivate
PERL5LIB="";
conda activate tormes-1.3.0
```
### sample sheet preparation
```bash
echo -e "Samples\tRead1\tRead2" >> samples_"$( date +"%Y-%m-%d" )".tab
for i in *R1_001.fastq.gz;do echo -e "${i/_S*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples_"$( date +"%Y-%m-%d" )".tab;done
```
### E.coli
```bash
tormes --metadata samples_"$( date +"%Y-%m-%d" )".tab --output ecoli_TORMES-"$( date +"%Y-%m-%d" )" --custom_genes_db ecoli_virulence --threads 32 --min_len 36 --genera Escherichia
```
### Salmonella
```bash
tormes --metadata samples.tab --output Salmonella_TORMES_2021 --threads 32 --min_len 36 --genera Salmonella
```
# Staphylococcus aureus analysis with bactopia
```bash
conda activate bactopia-new
for i in fastq/*fastq.gz;do mv $i ${i/_001};done
bactopia prepare fastq/ > fastqs.txt
bactopia --samples fastqs.txt -profile docker --datasets ./datasets --species "Staphylococcus aureus" --outdir Staphylococcus-aureus-out
bactopia --wf staphtyper --bactopia Staphylococcus-aureus-out -profile docker
bactopia --wf pangenome --bactopia Staphylococcus-aureus-out -profile docker
bactopia --wf mlst --bactopia Staphylococcus-aureus-out -profile docker
bactopia --wf abricate --bactopia Staphylococcus-aureus-out -profile docker
```
# fastANI analysis
```bash
conda activate assembly
for i in *fasta;do echo -e "`pwd`/$i" >> query_list;done
for i in *fasta;do echo -e "`pwd`/$i" >> reference_list;done
fastANI --ql query_list --rl reference_list --matrix -o fastani-out
sed -i "s/\/.*\///g" fastani-out.matrix
```

# Mycoplasma project
### kaken analysis
```bash
for i in A6_S5_L001_R1_001.fastq.gz; do kraken2 --db /home/vsingh/softwares/minikrake2_db/minikraken_8GB_20200312/ $i --threads 30 --output `basename ${i%_R1_001.fastq.gz}_kraken`;done
```
### Extracting reads according to taxanomy (544448)
```bash
git clone https://github.com/jenniferlu717/KrakenTools
```
```bash
for i in *fastq.gz;do python extract_kraken_reads.py -k ${i/L001*/L001_kraken} -s $i -t 544448 --include-children --report ${i/L001*/L001_R1_001_kraken2.report} -o ${i/R1/T544448_R1};done
```
### Assembly with spades
```bash
for i in *gz;do spades.py --only-assembler -t 50 -m 500 -k 21,33,55,77,99,121 -s $i -o ${i/_R1*/.ASM.out};done
```
### Tormes anslysis
```bash
conda deactivate
PERL5LIB="";
conda activate tormes-1.3.0
```
#### sample sheet preparation
```bash
echo -e "Samples\tRead1\tRead2" >> samples_"$( date +"%Y-%m-%d" )".tab
for i in *fasta;do echo -e "${i/.fasta*}\tGENOME\t`pwd`/$i" >> samples.tab;done
```

```bash
tormes --metadata samples_"$( date +"%Y-%m-%d" )".tab --output myco_TORMES-"$( date +"%Y-%m-%d" )" --threads 32
```



# Download SRA sequence read sample files from NCBI 
#### https://github.com/louiejtaylor/grabseqs
```bash
grabseqs sra ERR1777403
```
# Generating SNP density map chromosome wise from vcf file  

##convert vcf to simple 3 colums file, nano header SNP, Chr, Pos
### input vcf file "input.vcf"
```bash
more input.vcf | grep -v "#" | cut -f1,2,3 | awk '{ $3=$2; print $0}' | awk '{ $2=$1; print $0}' | awk '{ $1=$2__$3; print $0}' > input.txt
```

## R commands
```bash
library(CMplot)
PBR <-read.table("input.txt", header = TRUE)
CMplot(PBR,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300)
```

# snp density with ggplot2
```bash
snps<-read.table("SNPs.vcf",sep="\t",header=F,blank.lines.skip=TRUE,
                 comment.char = "#")
colnames(snps)<-c("chr","start","id","refallele","altallele","qual",
                  "filter","info","format")
summary(snps)
```
##### put the chromosomes in the good order: chr1, chr2, chr22, chrX
```bash
goodChrOrder <- c(paste("Supercontig_1.",c(1:8),sep=""),"MT_CBS_6936")
snps$chr <- factor(snps$chr,levels=goodChrOrder)
```

##### Plot the densities of snps in the bed file for each chr seperately
```bash
library(ggplot2)
snpDensity<-ggplot(snps) + 
geom_histogram(aes(x=start),binwidth=1e4) + # pick a binwidth that is not too small 
facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
ggtitle("Density of SNPs") +
xlab("Position in the genome") + 
ylab("SNP density") + 
theme_bw() # I prefer the black and white theme
```

##### save the plot to .pdf file
```bash
pdf("snp_density.pdf",800,1000)
print(snpDensity)
#dev.off()
```



