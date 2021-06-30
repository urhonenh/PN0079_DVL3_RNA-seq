#!/usr/bin/env Rscript

library(biomaRt)
library(dplyr)

# Ensembl gene IDs should be row names of the table.
ensemblgeneid_to_mgi_symbol <- function(ensemblgeneid_vector, prot_coding_only=FALSE) {
  
  mouse_ensembl <-  useMart("ensembl", dataset="mmusculus_gene_ensembl")
  
  gi <- getBM(attributes=c("chromosome_name", "start_position", "end_position", 
                           "strand", "ensembl_gene_id", "mgi_symbol", "gene_biotype"),
              filters="ensembl_gene_id", 
              values=as.character(ensemblgeneid_vector),
              mart=mouse_ensembl)
  
  # Select chromosomes.
  chr <- c(seq(1,22), "MT", "X", "Y")
  gi <- subset(gi, chromosome_name %in% chr)
  
  # Select only protein coding genes if desirable.
  if (prot_coding_only == TRUE) {
    
    gi <- gi[ gi$gene_biotype == 'protein_coding', ]
    
  }
    
  print(ensemblgeneid_vector[! ensemblgeneid_vector %in% gi$ensembl_gene_id])
  
  source("mousegenes_to_human.R")
  human_equival <- mousegenes_to_human(gi$mgi_symbol)
  gi <- merge(gi, human_equival, by.x="mgi_symbol", by.y="mgi_symbol", all.x = T)
  
  # Remove LRG ID rows.
  # is this needed? TEST THIS!
  # df[!grepl("^LRG", df$Name),]
  
  # Remove gene biotype column.
  # gi <- dplyr::select(gi, -c("gene_biotype"))
  
  return(gi)
}
