library("integration")
library("metagenomeSeq")
library("globaltest") 
library("SummarizedExperiment")
library("vegan") # For all the matrice

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
# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Check metadata with the names present in both datas
meta_r <- droplevels(meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ])

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

expr_norm <- norm_RNAseq(expr)
expr <- norm_RNAseq(expr_norm)
otus_table_i <- norm_otus(otus_table_i)

su <- apply(meta_r[, c("ID", "Exact_location", "AGE_SAMPLE", "diagTime", "SESCD_local", "SEX")], 2, is.na)
pos <- which(rowSums(su) == 1, arr.ind = TRUE)

meta_r$SEX <- as.factor(meta_r$SEX)

# Only for one gene
out <- lapply(seq_len(5000), function(x){
  gt(t(expr)[-pos, x]~Exact_location * IBD*ID*SEX + AGE_SAMPLE + diagTime, data = meta_r[-pos, ])
})

# Globaltest can be used in summarized experiment to asses if a whole SE is significant for a given relationship

meta_r$Exact_location <- as.factor(meta_r$Exact_location)
meta_r$Treatment <- as.factor(meta_r$Treatment)
meta_r$Transplant <- as.factor(meta_r$Transplant)
meta_r$Surgery <- as.factor(meta_r$Surgery)
phenoData <- AnnotatedDataFrame(droplevels(meta_r))
rownames(phenoData) <- colnames(expr)
rna <- ExpressionSet(expr, phenoData = phenoData)
otus <- as.matrix(otus_table_i)
colnames(otus) <- colnames(expr)
micro <- ExpressionSet(otus, phenoData = phenoData)

gt(SEX, rna)
gt(SEX, micro)
gt(Exact_location, rna[, !is.na(rna$Exact_location)])
gt(Exact_location, micro[, !is.na(micro$Exact_location)])
gt(diagTime, rna)
gt(diagTime, micro)
gt(AgeDiag, rna[, !is.na(rna$AgeDiag)])
gt(AgeDiag, micro[, !is.na(micro$AgeDiag)])
gt(Surgery, rna[, !is.na(rna$Surgery)])
gt(Surgery, micro[, !is.na(micro$Surgery)])
gt(Transplant, rna)
gt(Transplant, micro)
gt(Transplant:SEX, rna)
gt(Transplant:SEX, micro)
gt(Transplant:SEX:Surgery, rna)
##   p-value Statistic Expected Std.dev  #Cov
##1 1.34e-06      1.59    0.649  0.0888 37662
gt(Transplant:SEX:Surgery, micro)

gt(Transplant:SEX:Surgery:Exact_location, rna[, !is.na(rna$Exact_location)])
##    p-value Statistic Expected Std.dev  #Cov
## 1 2.55e-72      1.95    0.658  0.0419 37662
gt(Transplant:SEX:Surgery:Exact_location, micro[, !is.na(micro$Exact_location)])
##    p-value Statistic Expected Std.dev #Cov
## 1 4.08e-10      1.14    0.658  0.0458  529

saveRDS(rna, "rna_ES.RDS")
saveRDS(micro, "micro_ES.RDS")

all_RNA <- adonis(as.data.frame(t(expr))[-pos, ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[-pos, ], method = "euclidean", permutations = 5000)
all_16S <- adonis(as.data.frame(t(otus_table_i))[-pos, ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[-pos, ], method = "euclidean", by = "margin", permutations = 5000)

# CONTROLS
contr <- setdiff(which(meta_r$IBD == "CONTROL"), pos)

c_RNA <- adonis(as.data.frame(t(expr))[contr, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[contr, ], method = "euclidean")
c_16S <- adonis(as.data.frame(t(otus_table_i))[contr, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[contr, ], method = "euclidean", by = "margin")


# IBD
ibd <- setdiff(which(meta_r$IBD != "CONTROL"), pos)

i_RNA <- adonis(as.data.frame(t(expr))[ibd, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local, meta_r[ibd, ], method = "euclidean")
i_16S <- adonis(as.data.frame(t(otus_table_i))[ibd, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local, meta_r[ibd, ], method = "euclidean", by = "margin")

save(all_RNA, all_16S, c_RNA, c_16S, i_RNA, i_16S, file = "models.RData")
