cd <- setwd("..")

# Load the helper file
source("helper_functions.R")

intestinal <- "intestinal_16S"

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Load the input data
rna <- "intestinal_RNAseq"
expr <- read.delim(file.path(rna, "table.counts.results"), check.names = FALSE)
file_meta_r <- file.path(rna, "111217_metadata.csv")
meta_r <- read.table(file_meta_r, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = ";",
                     na.strings = c(NA, ""), header = TRUE, dec = c(",", "."))

setwd(cd)

# Summarize to genus
library("metagenomeSeq")

# Create the objects to summarize data
MR_i <- newMRexperiment(otus_table_i, 
                        # phenoData = AnnotatedDataFrame(meta_r),
                        featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i)))
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")

# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# Find the samples that we have microbiota and expression
int <- intersect(meta_r$Sample_Code_uDNA[!is.na(meta_r$Sample_Code_uDNA) & 
                                           !is.na(meta_r$`Sample Name_RNA`)],
                 meta_i$Sample_Code)

meta_i <- meta_i[meta_i$Sample_Code %in% int, ]
meta_r <- meta_r[meta_r$Sample_Code_uDNA %in% int, ]
meta_r <- meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Match the labels and order to append the id
meta_i <- meta_i[match(meta_r$Sample_Code_uDNA, meta_i$Sample_Code), ]
meta_r$`Sample Name_Code` <- gsub("([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1", rownames(meta_i))

colnames(genus_i) <- gsub("([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1", 
                          colnames(genus_i))

# Subset expression and outs
expr <- expr[, meta_r$`Sample Name_RNA`]
genus_i <- genus_i[, meta_r$`Sample Name_Code`]

# Subset if all the rows are 0 and if sd is 0
genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ] 
expr <- expr[apply(expr, 1, sd) != 0, ] 

abundance <- 0.005 # 0.5%

## All samples ####
# Filter by abundance at 0.5%
a <- prop.table(genus_i, 2)
b <- rowSums(a > abundance)

genus_i <- genus_i[b != 0, ]
# Correlate
p <- cor(t(genus_i), t(expr))
saveRDS(p, file = "correlations.RDS")

## All IBD ####
disease_i <- genus_i[, meta_r$IBD == "CD"]
disease_r <- expr[, meta_r$IBD == "CD"]

disease_i <- disease_i[apply(disease_i, 1, sd) != 0, ] 
disease_r <- disease_r[apply(disease_r, 1, sd) != 0, ] 

# Filter by abundance at 0.5%
a <- prop.table(disease_i, 2)
b <- rowSums(a > abundance)

disease_i <- disease_i[b != 0, ]

# Correlate
p <- cor(t(disease_i), t(disease_r))
saveRDS(p, file = "correlations_IBD.RDS")


## All Controls ####
disease_i <- genus_i[, meta_r$IBD == "CONTROL"]
disease_r <- expr[, meta_r$IBD == "CONTROL"]

disease_i <- disease_i[apply(disease_i, 1, sd) != 0, ] 
disease_r <- disease_r[apply(disease_r, 1, sd) != 0, ] 

# Filter by abundance at 0.5%
a <- prop.table(disease_i, 2)
b <- rowSums(a > abundance)

disease_i <- disease_i[b != 0, ]

# Correlate
p <- cor(t(disease_i), t(disease_r))
saveRDS(p, file = "correlations_C.RDS")
