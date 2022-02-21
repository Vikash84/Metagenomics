# This exercise is taken from the excellent Dada2 tutorial
# You can read that tutorial in depth here: https://benjjneb.github.io/dada2/tutorial.html
# This tutorial was written by: Benjamin Callahan 
# It is distributed under a CC BY 4.0 license

# For convienience, I've copied all the commands you need to download the packages
# here and have copied the relevant code chunks. 

# To get this tutorial to run on the lab computer or your own computer 
# you will need a copy of dada2 and phyloseq (as you know this is sometimes the most time
# consuming part of the whole process). 

# You will also need the example data to complete the tutorial. 
# The links to download that data are included here a few lines down. 


# We will run through this example in class today.
# There is no lab report associated with this tutorial.

# I don't expect you to perfectly understand all the commands that are run in this tutorial
# However, you should be able to understand WHY the steps are being completed
# and generally what the steps are doing.

# If you have questions, please consult the original tutorial for in-depth explainations 
# as you review this script after class.

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions

# Follow the promts to update packages if needed.

# The install will take a few minutes. 

# load the pacakge:
library(dada2); packageVersion("dada2")

# download the example files: http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip

# store a variable named path
path <- "~/Desktop/MiSeq_SOP/"
# set the path
setwd(path)


# take a look at the files in the folder:
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# visualize the quality profiles of forward reads
plotQualityProfile(fnFs[1:2])

# These look good. 
# We will truncate the forward reads at position 240 (trimming the last 10 nucleotides).
# Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs[1:2])

#  Based on these profiles, we will truncate the reverse reads at position 160 
# where the quality distribution crashes.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn error rates for forward and reverse reads
# These commands will take about 3 minutes each to run

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Apply the core sample infrence algorithm for F and R
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# merge the paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merged files
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Make a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Chimeras account for about 4% of the merged reads.


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# download the silva_nr_v132_train_set.fa.gz file from this link:
# https://zenodo.org/record/1172783#.XdwdeOdKjUI
# put it in the same directory as the sequence files

# Assign taxonomy
# this will take a pretty long time
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", multithread=TRUE)

# inspect the assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Evaluate accuaracy of DADA2 on the mock community of 20 known strains

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


# One nice thing about doing assingments with DADA2 is that it plays nicely with phyloseq
# Like dada2 this doesn't exisit in CRAN, so instalation is a bit tricky
# This can take a while, it's a big package.

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")


theme_set(theme_bw())

# Make a sample data table from file names
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


# construct a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


# now we can use some of the built in functions in phyloseq to look at our data
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")


# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# make a bar plot

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
