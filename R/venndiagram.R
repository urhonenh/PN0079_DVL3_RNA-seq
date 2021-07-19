#!/usr/bin/env Rscript
# A script for plotting Venn diagrams from RNA-seq DEG lists (created with DESeq2).

suppressPackageStartupMessages({
  library(VennDiagram)
  library(RColorBrewer)
})

args = commandArgs(trailingOnly = TRUE)
de_filename_1 <- args[1]
de_filename_2 <- args[2]
de_filename_3 <- args[3]
upregulated <- args[4]  # TRUE/FALSE

################

# Without last "/"
dedir <- "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/04_DE_analysis/coding_genes_only/cpm_1_filtered_genes"
outdir <- "/bmt-data/genomics/projects/dvl3_mouse_rna-seq/06_visualization"

if (dir.exists(paste0(outdir, "/venn_diagrams_gene_lists")) == FALSE) {
  dir.create(paste0(outdir, "/venn_diagrams_gene_lists"))  # Make a directory for the Venn diagram gene lists.
}

myCol <- brewer.pal(3, "Pastel1")

if (de_filename_3 != "-"){  # Three comparisons to one Venn diagram
  
  deg1 <- read.table(paste0(dedir, "/", de_filename_1), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  deg2 <- read.table(paste0(dedir, "/", de_filename_2), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  deg3 <- read.table(paste0(dedir, "/", de_filename_3), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  
  if (upregulated == TRUE){  # Up-regulated genes
    
    direction = "Up-regulated"
    mainpos = c(0.5, 0.62)  # 0.62 / 0.7
    g1 <- deg1[(deg1$log2FoldChange > 0), "mgi_symbol"]  # There were some missing symbols -> use Ensembl gene ids instead?
    g2 <- deg2[(deg2$log2FoldChange > 0), "mgi_symbol"]
    g3 <- deg3[(deg3$log2FoldChange > 0), "mgi_symbol"]
    
  } else {   # Down-regulated genes
    
    direction = "Down-regulated"
    mainpos = c(0.5, 0.62)  # 0.7
    g1 <- deg1[(deg1$log2FoldChange < 0), "mgi_symbol"]
    g2 <- deg2[(deg2$log2FoldChange < 0), "mgi_symbol"]
    g3 <- deg3[(deg3$log2FoldChange < 0), "mgi_symbol"]
  }
  
  g1 <- g1[!(is.na(g1) | g1=="")]
  g2 <- g2[!(is.na(g2) | g2=="")]
  g3 <- g3[!(is.na(g3) | g3=="")]
  
  de_filename_1 <- gsub("_padj_0.05_log2fc_1.5.txt", "", de_filename_1)
  de_filename_2 <- gsub("_padj_0.05_log2fc_1.5.txt", "", de_filename_2)
  de_filename_3 <- gsub("_padj_0.05_log2fc_1.5.txt", "", de_filename_3)
  
  de_name_1 <- gsub("_", " ", de_filename_1)
  de_name_2 <- gsub("_", " ", de_filename_2)
  de_name_3 <- gsub("_", " ", de_filename_3)
  
  venn.plot <- venn.diagram(x = list(g1, g2, g3),
                            main = paste0(direction, " DE genes (protein coding)"),
                            main.cex = 3,
                            main.pos = mainpos,
                            category.names = c(paste0(de_name_1, "\n(1)"), paste0(de_name_2, "\n(2)"), 
                                               paste0(de_name_3, "\n(3)")),
                            filename = NULL,
                            fill = myCol,
                            cex=3,
                            cat.cex = 2,
                            cat.dist = 0.05,
                            lwd = 2,
                            alpha = 0.6,
                            margin = 0.6,
                            # euler.d = TRUE,
                            # scaled = TRUE,
                            #sep.dist = 1,
                            rotation = 3,
                            ext.text = TRUE,)

  pdf(paste0(outdir, "/", de_filename_1, "+", de_filename_2, "+", de_filename_3, "_", direction, 
             "_coding_genes_only_venndiagram.pdf"), width=19, height=17)
  grid.newpage()
  grid.draw(venn.plot)
  dev.off()
  
  overlaps <- calculate.overlap(x = list(g1, g2, g3))
  names(overlaps) <- c("123", "12", "13", "23", "1", "2", "3")
  
  for (i in 1:length(overlaps)){
    
    o <- overlaps[[i]]
    write.table(o, paste0(outdir, "/venn_diagrams_gene_lists/", de_filename_1, "+", de_filename_2, "+",
                             de_filename_3, "_", direction, "_venn_overlaps_", names(overlaps)[[i]], "_coding_genes_only.txt"), 
                              sep = "\t", quote=F, col.names=F, row.names=F)
  }
  
} else {  # Two comparisons
  
  deg1 <- read.table(paste0(dedir, "/", de_filename_1), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  deg2 <- read.table(paste0(dedir, "/", de_filename_2), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  
  if (upregulated == TRUE){  # Up-regulated genes
    
    direction = "Up-regulated"
    g1 <- deg1[(deg1$log2FoldChange > 0), "mgi_symbol"]
    g2 <- deg2[(deg2$log2FoldChange > 0), "mgi_symbol"]
    
  } else {   # Down-regulated genes
    
    direction = "Down-regulated"
    g1 <- deg1[(deg1$log2FoldChange < 0), "mgi_symbol"]
    g2 <- deg2[(deg2$log2FoldChange < 0), "mgi_symbol"]
  }
  
  g1 <- g1[!(is.na(g1) | g1=="")]
  g2 <- g2[!(is.na(g2) | g2=="")]
  
  de_filename_1 <- gsub("_padj_0.05_log2fc_1.5.txt", "", de_filename_1)
  de_filename_2 <- gsub("_padj_0.05_log2fc_1.5.txt", "", de_filename_2)
  
  de_name_1 <- gsub("_", " ", de_filename_1)
  de_name_2 <- gsub("_", " ", de_filename_2)
  
  venn.plot <- venn.diagram(x = list(g1, g2),
                            main = paste0(direction, " DE genes (protein coding)"),
                            main.cex = 3,
                            main.pos = c(0.5, 0.65),
                            category.names = c(paste0(de_name_1, "\n(1)"), paste0(de_name_2, "\n(2)")),
                            filename = NULL,
                            fill = myCol[1:2],
                            cex=3,
                            cat.cex = 2,
                            cat.dist = 0.05,
                            #cat.pos = c(45, -45),
                            lwd = 2,
                            alpha = 0.6,
                            margin = 0.7)

  pdf(paste0(outdir, "/", de_filename_1, "+", de_filename_2, "_", direction, 
             "_coding_genes_only_venndiagram.pdf"), width=16, height=14)
  pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))
  grid.draw(venn.plot)
  dev.off()
  
  overlaps <- calculate.overlap(x = list(g1, g2))
  names(overlaps) <- c("12", "1", "2")
  
  for (i in 1:length(overlaps)){
    
    o <- overlaps[[i]]
    write.table(o, paste0(outdir, "/venn_diagrams_gene_lists/", de_filename_1, "+", de_filename_2, "_", direction, 
                "_venn_overlaps_", names(overlaps)[[i]], "_coding_genes_only.txt"), 
                sep = "\t", quote=F, col.names=F, row.names=F)
  }
}