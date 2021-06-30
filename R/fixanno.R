#!/usr/bin/env Rscript
# A function for fixing gene annotation columns and reordering the count table.
# Written by: Henna Urhonen
# E-mail: henna.urhonen@tuni.fi

library(dplyr)

# Input: A count table with Ensembl gene IDs as row names.
fixanno <- function(counttable) {
  
  # Remove ones from the strand column.
  counttable$strand <- gsub("-1", "-", counttable$strand)
  counttable$strand <- gsub("1", "+", counttable$strand)
  
  names(counttable)[names(counttable) == "Row.names"] <- "ensembl_gene_id"
  
  # Add "chr" to chromosome names if it is not included already.
  if (! "chr" %in% counttable$chromosome_name){
    counttable$chromosome_name <- paste0("chr", counttable$chromosome_name)
  }
  
  # Order chromosomes.
  chr_order <-c(paste0("chr", c((1:22),"MT","X","Y")))
  counttable$chromosome_name <- factor(counttable$chromosome_name, chr_order, ordered=TRUE)
  
  # Order table by given columns.
  counttable <- counttable[do.call(order, counttable[, c("chromosome_name", "start_position")]), ]
  
  # Reorder columns with dplyr select(). Use package name to avoid error.
  counttable <- counttable %>%
    dplyr::select("chromosome_name", "start_position", "end_position", 
                  "strand", "ensembl_gene_id", "mgi_symbol", "hgnc_symbol", "gene_biotype", everything())  # Note: mgi_symbol!
  
  return(counttable)
}