#!/usr/bin/env Rscript
# A script for running pairwise comparisons for DVL3 RNA-seq dataset with DESeq2.
# Written by: Henna Urhonen
# E-mail: henna.urhonen@tuni.fi

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

args = commandArgs(trailingOnly = TRUE)

# IMPORTANT! Remember to convert padj_threshold to numeric when ever it's used as a number!
padj_threshold <- as.numeric(args[1])  # Adjusted p-value threshold for filtering DESeq2 results.
log2fc_threshold <- as.numeric(args[2])   # Log2FoldChange threshold for filtering DESeq2 results.

strand_column=2  # 1=unstranded, 2=forward-stranded, 3=reverse-stranded (check!)

alignmentdir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/03_alignments/trimmed_bbduk"
outdir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/04_DE_analysis/"
gitrepo = "/home/hu425279/PN0079_DVL3_rna-seq/"

################

filepaths=list.files(path=alignmentdir, pattern="*ReadsPerGene.out.tab", full.names=TRUE)
samplenames=list.files(path=alignmentdir, pattern="*ReadsPerGene.out.tab", full.names=FALSE)

# Load a function for converting Ensembl gene IDs to gene symbols. NOTE! Mouse genome!
source ("ensemblgeneid_to_mgi_symbol.R")
source("fixanno.R")

# Remove parts of file names to make sample names match those in the colData file.
samplenames <- gsub('_ReadsPerGene.out.tab', '', samplenames)
samplenames <- gsub('_S[0-9][0-9]', '', samplenames)
samplenames <- gsub('_S[0-9]', '', samplenames)

# Set Ensembl gene IDs as row names and name the list elements according to the samples.
count_files = lapply(filepaths, function(x) read.table(x, header=F, skip=4, row.names=1))
names(count_files) <- samplenames

coldata <- read.table(paste(gitrepo, "data/coldata_deseq2.txt", sep=""), header=T, sep="\t", 
                      row.names=1, stringsAsFactors=T)
print(as.data.frame(coldata))
print(lapply(coldata, class))

# Add the raw read counts to a matrix.
mat <- matrix(data = NA, nrow=nrow(count_files[[1]]), ncol=nrow(coldata))
ids <- rownames(coldata)

for (i in 1:length(ids)) {
  sample <- ids[i]
  counts <- count_files[[sample]]
  counts <- counts[, strand_column] # Keep only the column with forward-stranded counts. Check the column!
  mat[,i] <- counts
}

colnames(mat) <- ids
rownames(mat) <- rownames(count_files[[1]])

# DESeq2 should also produce an error if these two vectors don't match.
if(all(rownames(coldata) != colnames(mat)) || !identical(rownames(coldata), colnames(mat)) || !all(rownames(coldata) %in% colnames(mat))) {
  stop("Row names of coldata don't match with column names of count matrix.")
}

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = coldata,
                              design = ~ Experiment_type + Radiation)  # Simple design for object construction. Remember to replace!

# ddsMat_LRT <- DESeq(dds, test="LRT", reduced =~ Time_point + Radiation + Experiment_type,
#                     full =~ Time_point + Radiation + Experiment_type + Radiation:Experiment_type)
# 
# design(dds) <- ~ Time_point + Experiment_type + Radiation + Experiment_type:Radiation

dds$group <- factor(paste0(dds$Experiment_type, "_", dds$Radiation))
design(dds) <- ~ group
dds$group <- relevel(dds$group, "Cell_line_0")
print(levels(dds$group))

write.table(model.matrix(object = design(dds), data = data.frame(colData(dds))),
            paste0(outdir, "deseq2_model_matrix_radiation+experiment_type.txt"), quote=F, sep="\t")

# Run DESeq2.
# NOTE! No pre-filtering of low count genes before DESeq2 for this QuantSeq dataset because lower counts are common. 
dds <- DESeq(dds)

saveRDS(dds, paste0(outdir, "dds_object.rds"))

resnames <- read.table(paste0(gitrepo, "data/de_comparisons.txt"), header=F, sep="\t", stringsAsFactors=F)

# Function for removing NAs.
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Pair-wise comparisons.
#for (i in 2:length(resnames)) {
for (j in 1:nrow(resnames)) {
  
  res <- results(dds, contrast=c("group", resnames[j,1], resnames[j,2]))
  
  comparison <- paste0(resnames[j,1], "_vs_", resnames[j,2])
  
  print(head(res))
  res <- as.data.frame(res)
  
  # Strip version numbers from Ensembl gene IDs.
  rownames(res) <- sub("\\..*.","\\.", rownames(res))
  rownames(res) <- gsub("\\.", "", rownames(res))
  
  # Remove NAs based on the adjusted p-value column.
  res <- completeFun(res, "padj")
  
  # Filter the table by log2FoldChange and adjusted p-value.
  res <- res[res$padj <= as.numeric(padj_threshold),]  # Remember as.numeric() !
  res <- res[abs(res$log2FoldChange) >= as.numeric(log2fc_threshold),]
  
  # Convert Ensembl gene IDs to gene symbols.
  # A corresponding gene symbol may not be found for some Ensembl gene IDs,
  # but those rows are still kept in the final results.
  if (nrow(res) != 0) {
    
    geneinfo <- ensemblgeneid_to_mgi_symbol(rownames(res))
    res <- merge(res, geneinfo, by.x="row.names", by.y="ensembl_gene_id", all.x = T)
    
    print(paste0("Number of genes before removing duplicates: ", nrow(res)))
    
    res <- res[!duplicated(res[,c("mgi_symbol", "log2FoldChange", "pvalue", "padj")]),]
    print(paste0("Number of genes after removing duplicates: ", nrow(res)))
  }
  
  res <- res[order(res$padj),]
  
  # Rearrange the columns and rename chromosome and strand values.
  if (nrow(res) != 0) {
    res <- fixanno(res)
  }
  
  # Possibility to check that the Ensembl gene ID -> HGNC symbol was OK.
  print(paste0(res[1,"ensembl_gene_id"], "=", res[1,"mgi_symbol"]))
  
  write.table(res, paste0(outdir, comparison, "_padj_", padj_threshold, 
                          "_log2fc_", log2fc_threshold, ".txt"), sep="\t", row.names=FALSE, quote=FALSE)
}