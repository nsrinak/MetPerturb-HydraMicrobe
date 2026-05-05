# DADA2 pipeline 

library(dada2)

# Due to high number of the samples. The analysis will follow DADA2 implementation to the big data via 
# this link: https://benjjneb.github.io/dada2/bigdata_paired.html
# The forward files and reverse files were expected to be in the different directory

pathF <- "/Forward"
pathR <- "/Reverse"

filtpathF <- file.path(pathF, "filtered") # to create sub-directory for files afte processing 
filtpathR <- file.path(pathR, "filtered")

fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(225,225), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

filtpathF <- "/Forward/filtered"
filtpathR <- "/Reverse/filtered"

# After filtering the quality will be checked to make sure good quality, again with FASTQC and MULTIQC

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Assumes filename = samplename_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab.before <- makeSequenceTable(mergers)

saveRDS(seqtab, "seqtab.rds")

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
print(dim(seqtab)) # dim =
print(sum(seqtab)/sum(seqtab.before)) # ratio =

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
tax_spec <- addSpecies(tax, "xxx")

# Write to disk
saveRDS(seqtab, "seqtab_final.rds")
saveRDS(tax, "tax_final.rds")