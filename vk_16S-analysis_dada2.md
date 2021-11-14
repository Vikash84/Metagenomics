
# 16S tutorial
#### Vikash Singh
#### November 14, 2021
#### reference: https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2019.10)
#### reference: https://docs.qiime2.org/2019.10/tutorials/

# 1 data preparation
# 1.0.1 copy/transfer all the fastq.gz files to fastq_files dir
```bash
mkdir fastq_files
```
# 1.0.2 make manifest file
######## job_manifest.tsv tab-separated values (TSV) file
######## ref: https://github.com/qiime2/docs/blob/master/source/tutorials/importing.rst#fastq-manifest-formats
```bash
sample-id   forward-absolute-filepath   reverse-absolute-filepath
Con5090Ileum    $PWD/fqgz/1183-1_S1_L001_R1_001.fastq.gz    $PWD/fqgz/1183-1_S1_L001_R2_001.fastq.gz
Con8124Ileum    $PWD/fqgz/1183-19_S19_L001_R1_001.fastq.gz  $PWD/fqgz/1183-19_S19_L001_R2_001.fastq.gz
...
```
