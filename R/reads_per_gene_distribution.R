#!/usr/bin/env Rscript
# A script for plotting distribution of reads per gene in the dataset.
# Written by: Henna Urhonen
# E-mail: henna.urhonen@tuni.fi

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

dedir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/04_DE_analysis/"
outdir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/06_visualization/"
gitrepo = "/home/hu425279/PN0079_DVL3_rna-seq/"

###############

dds <- readRDS(paste0(dedir, "dds_object.rds"))

raw_counts <- counts(dds)

raw_rowsums <- as.data.frame(rowSums(raw_counts))
log2_rowsums <- log2(raw_rowsums + 1)  # Use pseudocount 1

print(min(log2_rowsums))
print(max(log2_rowsums))

# Histogram of log2 transformed values for the whole dataset.
p <- ggplot(log2_rowsums, aes(x=log2_rowsums[,1])) + 
  geom_histogram(binwidth=1, colour="black", fill="white") + ggtitle("Total log2 read counts per gene in the DVL3 dataset") + 
  xlab("log2 read count") + ylab("Number of genes") + theme_bw()

pdf(paste0(outdir, "log2_reads_per_gene.pdf"))
p
dev.off()

# Histogram of raw counts for genes with smaller total counts.
library("ggforce")

low_count_rowsums <- as.data.frame(raw_rowsums[which(raw_rowsums[,1] <= 30),])

yRange = c(0,3000)

p2 <- ggplot(low_count_rowsums, aes(x=low_count_rowsums[,1])) + 
  geom_histogram(binwidth=1, colour="black", fill="white") + ggtitle("Low count genes with max. 30 reads in total") +
  xlab("Raw read count") + ylab("Number of genes") + theme_bw() + facet_zoom(ylim = c(0, 3000))  # Zoom to this range in an additional figure

pdf(paste0(outdir, "raw_reads_per_low_count_gene.pdf"))
p2
dev.off()