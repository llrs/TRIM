library("biomaRt")
ensembl  <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters <- listFilters(ensembl)
filters[grep("chr", filters$name), ]

attributes <- listAttributes(ensembl)
attributes[grep("gene_id", attributes$name), ]

bmY <- getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version', "hgnc_symbol"), 
      filters = 'chromosome_name', 
      values = "Y", 
      mart = ensembl)


bmX <- getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version', "hgnc_symbol"), 
            filters = 'chromosome_name', 
            values = "X", 
            mart = ensembl)

