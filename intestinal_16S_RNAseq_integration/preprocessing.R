cd <- setwd("..")
library("integration")
library("metagenomeSeq")

intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Load the input data
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(cd)

# Correct the swapped samples and match metadata
expr <- norm_expr_colnames(expr)
# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]



# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

ctrls_expr <- expr[, meta_r$IBD == "CONTROL"]
ctrls_otus_table_i <- otus_table_i[, meta_r$IBD == "CONTROL"]

IBD_expr <- expr[, meta_r$IBD != "CONTROL"]
IBD_otus_table_i <- otus_table_i[, meta_r$IBD != "CONTROL"]

expr_norm <- norm_RNAseq(expr)
ctrls_expr_norm <- norm_RNAseq(ctrls_expr)
IBD_expr_norm <- norm_RNAseq(IBD_expr)

# Filter expression
expr <- filter_RNAseq(expr_norm)
ctrls_expr <- filter_RNAseq(ctrls_expr_norm)
IBD_expr <- filter_RNAseq(IBD_expr_norm)

otus_table_i <- norm_otus(otus_table_i, otus_tax_i)
ctrls_otus_table_i <- norm_otus(ctrls_otus_table_i, otus_tax_i)
IBD_otus_table_i <- norm_otus(IBD_otus_table_i, otus_tax_i)

# Save
saveRDS(otus_table_i, "otus_table.RDS")
saveRDS(ctrls_otus_table_i, "ctrls_otus_table.RDS")
saveRDS(IBD_otus_table_i, "IBD_otus_table.RDS")
saveRDS(otus_tax_i, "otus_tax.RDS")
saveRDS(expr, "expr.RDS")
saveRDS(ctrls_expr, "ctrls_expr.RDS")
saveRDS(IBD_expr, "IBD_expr.RDS")
saveRDS(droplevels(meta_r), "meta.RDS")
