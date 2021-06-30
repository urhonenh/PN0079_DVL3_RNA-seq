#!/usr/bin/env Rscript
# A function for converting mouse gene symbols to human gene symbol equivalents.
# Written by: Henna Urhonen
# E-mail: henna.urhonen@tuni.fi

library(biomaRt)

# Basic function to convert mouse to human gene names
mousegenes_to_human <- function(mouse_genesyms){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  conv = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values=mouse_genesyms, 
                   mart=mouse, 
                   attributesL=c("hgnc_symbol"), 
                   martL=human, 
                   uniqueRows=T)
  
  #humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  colnames(conv) <- c("mgi_symbol", "hgnc_symbol")
  print(head(conv))
  return(conv)
}