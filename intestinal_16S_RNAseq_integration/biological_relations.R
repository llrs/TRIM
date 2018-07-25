# Test if the loadings of the genes are from a relevant pathways.
cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"


epithelium <- read.csv("epithelium.csv")
epithelium <- epithelium$Epithelium

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

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- colnames(expr)[rev(position)]
colnames(expr) <- toupper(colnames(expr))
#To match metadata
colnames(expr) <- gsub("16-TM29", "16-TM30", colnames(expr)) 

# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]
meta_r <- droplevels(meta_r)

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Translate into Entrez
epithelium <- epitheliumE(epithelium)

# Do the biological analysis of all the samples
all_bootstrap <- readRDS("bootstrap.RDS")
all_sgcca <- readRDS("sgcca.RDS")
write.csv(integration::weights(all_sgcca), file = "RNAseq_weight_all.csv", na = "", row.names = FALSE)
write.csv(integration::weights_otus(all_sgcca, otus_tax_i), file = "16S_weight_all.csv", na = "", row.names = FALSE)
biological_relationships(all_sgcca, all_bootstrap, "all", otus_tax_i, epithelium, today)

# Controls
controls_bootstrap <- readRDS("bootstrap_Controls.RDS")
controls_sgcca <- readRDS("Controls.RDS")
write.csv(integration::weights(controls_sgcca), file = "RNAseq_weight_controls.csv", na = "", row.names = FALSE)
write.csv(weights_otus(controls_sgcca, otus_tax_i), file = "16S_weight_controls.csv", na = "", row.names = FALSE)
biological_relationships(controls_sgcca, controls_bootstrap, "controls", otus_tax_i, epithelium, today)

# IBD
IBD_bootstrap <- readRDS("bootstrap_IBD.RDS")
IBD_sgcca <- readRDS("IBD.RDS")
write.csv(integration::weights(IBD_sgcca), file = "RNAseq_weight_IBD.csv", na = "", row.names = FALSE)
write.csv(weights_otus(IBD_sgcca, otus_tax_i), file = "16S_weight_IBD.csv", na = "", row.names = FALSE)
biological_relationships(IBD_sgcca, IBD_bootstrap, "IBD", otus_tax_i, epithelium, today)
