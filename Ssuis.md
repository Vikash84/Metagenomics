# H.parasuis Analysis

### sample sheet preparation
```bash
for i in *R1_001.fastq.gz;do echo -e "${i/_BTIS*}\t`pwd`/$i\t`pwd`/${i/R1/R2}" >> samples.tab;done
```
### Downloading refence sequenece from NCBI

### Running pipeline
```bash
nullarbor.pl --cpus 50 --run --trim --mlst hparasuis --treebuilder iqtree_slow --taxoner centrifuge --name Hparasuis_analysis --ref ref.fasta --input samples.tab --outdir results
```

### Serotyping
```bash
blastn -task blastn-short -query Hparasuis_sero_primers.fasta -db HPS.all_merged.fasta -out HPS.all_merged.blastn -outfmt "6 qseqid qlen qcovs stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs 5 -num_threads 50
```



# SSuis Analysis

```bash
nullarbor.pl --cpus 50 --run --trim --mlst ssuis --treebuilder iqtree_slow --taxoner centrifuge --name ssuis_analysis --ref Ssuis_ref.fa --input samples.tab --outdir results
```
## Ssuis serotyping
```bash
./Ssuis_serotypingPipeline.pl --fastq_directory /home/vsingh/vdl/Mor_Project_147/serotyping_github/SsuisSerotyping_pipeline/SsuisSerotyping_pipeline/data --forward _R1_001 --reverse _R2_001 --ends pe
```
