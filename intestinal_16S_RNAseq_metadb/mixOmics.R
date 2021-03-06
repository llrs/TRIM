library("ggforce")
library("mixOmics")
library("integration")
library("fgsea")
cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

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
meta_r <- droplevels(meta_r)

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes=TRUE, log = TRUE)

# Filter expression
expr <- filter_RNAseq(expr_norm)

# Normalize OTUS
otus_table_i <- norm_otus(otus_tax_i, otus_table_i)

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r
keepCol <- sapply(metadb, is.factor)
nam <- c(
  "Exact_location", # Segment of the sample
  # superseeded by SESCD 
  # "Active_area", # Health stage of the sample
  # "IBD", # Disease or control
  "AGE_SAMPLE", # Age
  # "diagTime", # Time with disease
  # Not really needed induced by diagTime and age sample
  "AgeDiag", # Age at which the disease was diagnositcated
  "Transplant", # Stage of the treatment
  "ID", # Patient
  # "SESCD_local", # Clinical score
  "Treatment", # Further complications
  "Surgery", # Up to surgery?
  "SEX" # Male/female
) 

metadb <- model_RGCCA(meta_r, nam) 
# Set metadb with a sigle variable with several options
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i))

lapply(A, dim)

design = matrix(0.1, ncol = length(A), nrow = length(A), 
                dimnames = list(names(A), names(A)))
diag(design) = 0


# Use mixomics to interact with the data

Y <- meta_r$Exact_location
Y[is.na(Y)] <- "ILEUM"
A <- lapply(A, function(x){
  rownames(x) <- meta_r$Seq_code_uDNA
  x
})

Ab <- lapply(A, function(x){
  x[meta_r$Seq_code_uDNA[meta_r$IBD == "CD"], ]
})
# Not possible by ID
Y <- Y[meta_r$IBD == "CD"]

sgccda.res = block.splsda(X = A, Y = Y, ncomp = 5, 
                          design = design)
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10)
plot(perf.diablo) 
