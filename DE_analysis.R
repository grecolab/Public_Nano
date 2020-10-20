rm(list=ls())

library(dplyr)
library(tibble)
library(readr)
library(stats)
library(ggplot2)
library(Rsubread)
library(NOISeq)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(grid)
library(gridExtra)
library(Rsamtools)
library(ggridges)
library(sva)


mainDir <- “/path/to/files/"
#batch <- "batch2/"
dataset <- "GSEXXXX” #type in the GEO dataset ID
setwd(paste0(mainDir, dataset))

source(“get_and_filter_counts_functions.R”)

ann_human <- "/home/ESPRESSO/antonio/rna-seq_tools/gencode.v35.annotation.gtf"
ann_mouse <- "/home/ESPRESSO/antonio/rna-seq_tools/gencode.vM25.annotation.gtf"
#gtf <- rtracklayer::import("/home/MOKA/antonio/RNA-SeqTools/indexes/gencode.v30lift37.basic.annotation.2.gtf")

phenodata <- read_delim(paste0(mainDir, "SEQ_meta/", "GSEXXXX_metadata.txt"), col_names = TRUE, delim = "\t")

sratable <- read_delim(paste0(mainDir, dataset, "/filereport_read_run_PRJNAXXXX_tsv.txt"), col_names = TRUE, delim = "\t")

#phenodata$subject_status <- gsub(phenodata$subject_status, pattern = " ", replacement = "_", fixed = TRUE)
table(phenodata$group)


#phenodata$time_point[which(is.na(phenodata$time_point))] <- "untreated"

bam.files <- list.files(path = paste0(mainDir, dataset), pattern = "_unique_sorted.bam")
sample_name <- sapply(strsplit(bam.files, "_", fixed=TRUE), "[", 1)
bam.files <- sapply(bam.files, function(x) paste0(mainDir, dataset, "/", x))

#bam_sublist <- split(bam.files, ceiling(seq_along(bam.files)/5))

### !!! Set the library layout !!! ###

total_raw_counts <- get_raw_counts(bam.files, annotation = ann_human, paired.end = TRUE, ncores = 45)

colnames(total_raw_counts[[1]]) <- sample_name
colnames(total_raw_counts[[1]]) <- sratable$sample_alias[match(colnames(total_raw_counts[[1]]), sratable$run_accession)]


write.table(total_raw_counts[[1]], file = paste0("raw_counts_matrix_", dataset, ".txt"), quote = FALSE, sep = "\t", col=NA)


get_gene_lengths <- total_raw_counts[[2]]

### Conditions ###


group <- phenodata$group
conditions <- factor(group)

### Low counts filtering ###

filtered_data <- filter_low_counts(counts.matrix = total_raw_counts[[1]], conditions = conditions, method = "proportion", depth=as.numeric(apply(total_raw_counts[[1]], 2, sum)), normalized = FALSE, p.adj = "fdr")


write.table(filtered_data, file = paste0("filtered_counts_matrix_", dataset, ".txt"), quote = FALSE, sep = "\t", col=NA)


### DESeq2 ###

myfactors <- phenodata


#myfactors$phenodata.tissue_type <- factor(myfactors$phenodata.tissue_type, levels = c("psoriatic_skin", "Uninvolved_skin", "normal_skin"))
rownames(myfactors) <- colnames(filtered_data)

ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = filtered_data,
                                 colData = myfactors,
                                 design = ~group)

dds <- DESeq2::DESeq(ddsMat)
DESeq2::counts(dds, normalized=TRUE)
write.table(DESeq2::counts(dds, normalized=TRUE), file = paste0("normalized_counts_matrix_", dataset, ".txt"), quote = FALSE, sep = "\t", row.names = TRUE, col=NA)
#resultsNames(dds)
res1 <- results(dds, contrast = c("group","Li_doped_Ag2S_320", "Control"), pAdjustMethod = "BH")
write.table(res1, file = paste0("DEG_results_DESeq2_", dataset, "_", “condition_1_vs_Control", ".txt"), row.names = TRUE, sep = "\t", quote = FALSE)

res1_adj <- res1[which(res1$padj<=0.05),]



write.table(res1_adj, file = paste0("DEG_results_DESeq2_", dataset, "_", “condition_1_vs_Control", "_filtered_fdr.txt"), row.names = TRUE, sep = "\t", quote = FALSE)

