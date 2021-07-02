#!/usr/bin/env Rscript
# A script for making PCA plots with DESeq2.
# Henna Urhonen
# E-mail: henna.urhonen@tuni.fi

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

alignmentdir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/03_alignments/trimmed_bbduk"
outdir = "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/06_visualization/"
gitrepo = "/home/hu425279/PN0079_DVL3_rna-seq/"

###############

# OR:
# dds <- readRDS(paste0(outdir, "plotpca_deseq2dataset_object.rds"))

filepaths=list.files(path=alignmentdir, pattern="*ReadsPerGene.out.tab", full.names=TRUE)
samplenames=list.files(path=alignmentdir, pattern="*ReadsPerGene.out.tab", full.names=FALSE)

# Remove parts of file names to make sample names match those in the colData file.
samplenames <- gsub('_ReadsPerGene.out.tab', '', samplenames)
samplenames <- gsub('_S[0-9][0-9]', '', samplenames)
samplenames <- gsub('_S[0-9]', '', samplenames)

# Set Ensembl gene IDs as row names and name the list elements according to the samples.
count_files = lapply(filepaths, function(x) read.table(x, header=F, skip=4, row.names=1))
names(count_files) <- samplenames

coldata <- read.table(paste(gitrepo, "data/coldata_deseq2.txt", sep=""),
                      header=T, sep="\t", row.names=1, stringsAsFactors=T)
print(as.data.frame(coldata))
print(lapply(coldata, class))  # Check that all colData columns are factors.

# Add the raw read counts to a matrix.
mat <- matrix(data = NA, nrow=nrow(count_files[[1]]), ncol=nrow(coldata))
ids <- rownames(coldata)

if(all(ids %in% samplenames)==FALSE){
  stop("Warning: Not all row names of colData are included in sample names")
}

for (i in 1:length(ids)) {
  sample <- ids[i]
  counts <- count_files[[sample]]
  counts <- counts[,2] # Keep only the column with forward-stranded counts. Check the column!
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

dds$group <- factor(paste0(dds$Experiment_type, "_", dds$Radiation))
design(dds) <- ~ group
dds$group <- relevel(dds$group, "Cell_line_0")

saveRDS(dds, paste0(outdir, "plotpca_deseq2dataset_obj_group_design.rds"))
#design(dds) <- ~ Batch + Culture + Treatment + Culture:Treatment
#dds$Treatment <- relevel(dds$Treatment, "etoh")

# I tested both blind=FALSE and blind=TRUE and there were no differences in the PCA plots.
vsd <- vst(dds, blind=FALSE)

### PCA plots.
# Intgroup (design variables) need to be specified.

pdf(paste0(outdir, "pca_test_experiment_type.pdf"))
plotPCA(vsd, "Experiment_type")
dev.off()
 
pcaData_rad <- plotPCA(vsd, intgroup="Radiation", returnData=TRUE)
percentVar <- round(100 * attr(pcaData_rad, "percentVar"))
options(repr.plot.width = 14, repr.plot.height = 8)
pcaplot_rad <- ggplot(pcaData_rad, aes(PC1, PC2, color=Radiation)) +
  #scale_color_brewer(palette="Blues") +
  scale_color_grey(start=0.8, end=0.2) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

pcaplot_rad <- pcaplot_rad + guides(color=guide_legend(title="Radiation (Gy)"))

pdf(paste0(outdir, "pca_test_radiation.pdf"))
pcaplot_rad
dev.off()

# pdf(paste0(outdir, "pca_test_cell_type.pdf"))
# plotPCA(vsd, "Cell_type")
# dev.off()

# pdf(paste0(outdir, "pca_time_point.pdf"))
# plotPCA(vsd, intgroup=c("Time_point"))
# dev.off()

# pdf(paste0(outdir, "pca_test.pdf"))
# plotPCA(vsd, intgroup=c("Cell_type", "Experiment_type", "Radiation"))
# dev.off()

require("ggrepel")
set.seed(42)

pcaData <- plotPCA(vsd, intgroup=c("Radiation", "Experiment_type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
options(repr.plot.width = 14, repr.plot.height = 8)
pcaplot <- ggplot(pcaData, aes(PC1, PC2, shape=Radiation, color=Experiment_type)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) # +coord_fixed()

pcaplot <- pcaplot + geom_text_repel(aes(label=rownames(colData(dds))), colour = "black", size=2)
#pcaplot + geom_text_repel(aes(label=Radiation), colour = "black", size=3)
pdf(paste0(outdir, "pca_radiation_and_experiment_type.pdf"))
pcaplot
dev.off()