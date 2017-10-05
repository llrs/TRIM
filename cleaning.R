otus_table_i <- read.csv("intestinal_16S/OTUs-Table-new-biopsies.csv", stringsAsFactors = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]

#' Clean and prepare the data from IMNGS
#' 
#' Divides the taxonomy into a new matrix for each otu
#' 
#' @param taxonomy Last column of files a string ; separated with domain, 
#' phylum, vlass, order, family, genus and species.
#' @param otus The name of the rows
#' 
#' @return 
#' A matrix with the taxonomic information ready for the package phylo
taxonomy <- function(taxonomy, otus){
  taxonomy <- sapply(taxonomy, strsplit, split = ";")
  names(taxonomy) <- otus
  otus_tax <- t(sapply(taxonomy, '[', seq(max(sapply(taxonomy, length)))))
  colnames(otus_tax) <- c("Domain", "Phylum", "Class", "Order", 
                          "Family", "Genus", "Species")
  # Remove spaces
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "\\s", replacement = "")
  otus_tax[otus_tax == ""] <- NA # Remove empty cells
  otus_tax
}

otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Read for the stools
otus_table_s <- read.delim("stools_16S/OTUs-Table-refined-stools.tab", stringsAsFactors = FALSE)
tax_s <- otus_table_s[, ncol(otus_table_s)]
otus_tax_s <- taxonomy(tax_s, rownames(otus_table_s))
