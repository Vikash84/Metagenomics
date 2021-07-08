# H.parasuis Analysis

### sample sheet preparation
```bash
for i in *R1_001.fastq.gz;do echo -e "${i/_BTIS*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples.tab;done
```
### Downloading refence sequenece from NCBI

### Running pipeline
```bash
PERL5LIB="";nullarbor.pl --cpus 50 --run --trim --mlst hparasuis --treebuilder iqtree_slow --taxoner centrifuge --name Hparasuis_analysis --ref ref.fasta --input samples.tab --outdir results
```

### Serotyping
```bash
blastn -task blastn-short -query Hparasuis_sero_primers.fasta -db HPS.all_merged.fasta -out HPS.all_merged.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50

sort -k3,3nr -k15,15nr HPS.all_merged.blastn | awk -F"\t" '!seen[$1]++' | less
```
### vtaA profile
```bash
blastn -task blastn-short -query vtaA_profile.fasta -db scaffolds.fasta -out vtaA_profile.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50

sort -k3,3nr -k15,15nr vtaA_profile.blastn | awk -F"\t" '!seen[$1]++' | less
```



# SSuis Analysis
```bash
PERL5LIB="";nullarbor.pl --cpus 50 --run --trim --mlst ssuis --treebuilder iqtree_slow --taxoner centrifuge --name ssuis_analysis --ref Ssuis_ref.fa --input samples.tab --outdir results
```
### Ssuis serotyping
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
conda deactivate base
PERL5LIB="";
conda activate tormes-1.3.0
```
### sample sheet preparation
```bash
echo -e "Samples\tRead1\tRead2" >> samples.tab
for i in *R1_001.fastq.gz;do echo -e "${i/_BTIS*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples.tab;done
```
### E.coli
```bash
tormes --metadata samples.tab --output ecoli_TORMES-070821 --custom_genes_db ecoli_virulence --threads 32 --genera Escherichia
```
### Salmonella
```bash
tormes --metadata samples.tab --output Salmonella_TORMES_2021 --threads 32 --genera Salmonella
```

# fastANI analysis
```bash
conda activate assembly
for i in *fasta;do echo -e "`pwd`/$i" >> query_list;done
for i in *fasta;do echo -e "`pwd`/$i" >> reference_list;done
fastANI --ql query_list --rl reference_list -o fastani-out
```


