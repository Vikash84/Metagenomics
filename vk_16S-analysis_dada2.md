# 16S tutorial
###### Vikash Singh
###### November 14, 2021
###### reference: https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2019.10)
###### reference: https://docs.qiime2.org/2019.10/tutorials/

# 1 data preparation
## 1.0.1 copy/transfer all the fastq.gz files to fastq_files dir
```bash
mkdir fastq_files
```
## 1.0.2 make manifest file
###### job_manifest.tsv tab-separated values (TSV) file
###### ref: https://github.com/qiime2/docs/blob/master/source/tutorials/importing.rst#fastq-manifest-formats
```bash
sample-id   forward-absolute-filepath   reverse-absolute-filepath
Con5090Ileum    $PWD/fqgz/1183-1_S1_L001_R1_001.fastq.gz    $PWD/fqgz/1183-1_S1_L001_R2_001.fastq.gz
Con8124Ileum    $PWD/fqgz/1183-19_S19_L001_R1_001.fastq.gz  $PWD/fqgz/1183-19_S19_L001_R2_001.fastq.gz
...
```
## 1.0.3 make metadata file
###### for whole job
###### job_meta.tsv tab-separated values (TSV) file
###### ref: https://docs.qiime2.org/2019.10/tutorials/metadata/
```bash
sample-id   GroupID treatment-group E.coliChallenge Sex Euth    PigID   Sourceofsample  Datetaken   NGS-SampleNo
#q2:types   categorical categorical categorical categorical categorical categorical categorical categorical categorical
Con5090Ileum    C   Control NO  M   1   5090    Ileum   1/7/2019    1
Con8124Ileum    C   Control NO  F   2   8124    Ileum   5/7/2019    19
Con8141Ileum    C   Control NO  M   2   8141    Ileum   5/7/2019    20
...
```

###### for categorical variables
###### e.g.GroupID tab-separated values (TSV) file
```bash
sample-id   GroupID
#q2:types   categorical
C   C
CR  CR
CR-EC   CR-EC
EC  EC
```
## 1.0.4 prepare/download database for taxonomy assignment
We use SILVA database **silva_132_99_V4/silva-132-99-515-806-nb-classifier.qza**

# 2 QIIME2 anlysis steps
## 2.1 data importing
###### reads type" PairedEndFastqManifestPhred33V2
###### code in qiime2:
```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path job_manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
the above code produces the **paired-end-demux.qza file

## 2.2 primer trimming with cutadapt
###### ref https://github.com/SchlossLab/MiSeq_WetLab_SOP

###### DNA was amplified by using the 515f/806r primer set:
###### Forward V4: GTGCCAGCMGCCGCGGTAA
###### Reverse V4: GGACTACHVGGGTWTCTAAT

code in qiime2:
```bash
qiime cutadapt trim-paired \
   --i-demultiplexed-sequences paired-end-demux.qza \
   --p-cores 20 \
   --p-front-f GTGCCAGCMGCCGCGGTAA \
   --p-front-r GGACTACHVGGGTWTCTAAT \
   --o-trimmed-sequences pe_reads_cutadapt_trimmed.qza \
   --verbose \
   &> primer_trim.log 
```
--p-cores is number of cores. the above code produces the **pe_reads_cutadapt_trimmed.qza file.
## 2.3 Summarize trimmed FASTQs
###### Check quality plots and sequence length
###### code in qiime2:
```bash
qiime demux summarize \
    --i-data pe_reads_cutadapt_trimmed.qza \
    --o-visualization pe_reads_cutadapt_trimmed.qzv
```
##### To view and check the pe_reads_cutadapt_trimmed.qzv, you can
<!--- qiime tools view pe_reads_cutadapt_trimmed.qzv , or
upload this file to view.qiime2.org, or
unzip the qzv and open the data/index.html.--->

Based on the plots you see in qzv, decide values would you choose for --p-trunc-len and --p-trim-left in DADA2 denoising ref: https://docs.qiime2.org/2019.10/tutorials/atacama-soils/

