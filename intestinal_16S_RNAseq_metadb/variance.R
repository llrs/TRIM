cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

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
file_meta_r <- file.path(rna, "metadata_13032018.csv")
meta_r <- read.table(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, sep = ";",
  na.strings = c(NA, ""), header = TRUE, dec = c(",", ".")
)

setwd(cd)

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- colnames(expr)[rev(position)]

# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# Find the samples that we have microbiota and expression in the metadata
int <- intersect(meta_r$Sample_Code_uDNA, meta_i$Sample_Code)

# Subset metadata
meta_i <- meta_i[meta_i$Sample_Code %in% int, ]
meta_r <- meta_r[meta_r$Sample_Code_uDNA %in% int, ]
# Duplicate label! the one we don't know if it is colon or ileum

# Match the labels and order to append the id
meta_i <- meta_i[match(meta_r$Sample_Code_uDNA, meta_i$Sample_Code), ]

# Add sample name
pattern <- "([0-9]{2,3}\\.B[0-9]+)\\..+"

meta_r$`Sample Name_Code` <- gsub(pattern, "\\1", rownames(meta_i))
colnames(otus_table_i) <- gsub(pattern, "\\1", colnames(otus_table_i))

# Intersect between sequencing and metadata
int_16S <- intersect(meta_r$`Sample Name_Code`, colnames(otus_table_i))
int_RNAseq <- intersect(meta_r$`Sample Name_RNA`, colnames(expr))

# Check that metadata for the samples match
int2 <- intersect(
  meta_r[
    meta_r$`Sample Name_RNA` %in% int_RNAseq,
    "Sample_Code_uDNA"
    ],
  meta_i[int_16S, "Sample_Code"]
)
# Subset metadata
meta_i <- meta_i[meta_i$Sample_Code %in% int2, ]
meta_r <- meta_r[meta_r$Sample_Code_uDNA %in% int2, ]
# Still carring the duplicate label!

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$`Sample Name_Code`]

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes=TRUE, log = TRUE)

# Filter expression
expr <- norm_RNAseq(expr_norm)

# Normalize OTUS
library("metagenomeSeq")
MR_i <- newMRexperiment(
  otus_table_i, 
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i[rownames(otus_table_i), ]))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus_table_i <- MRcounts(MR_i, norm = TRUE, log = TRUE)

# Subset if all the rows are 0 and if sd is 0
otus_table_i <- otus_table_i[apply(otus_table_i, 1, sd) != 0, ]

library("globaltest") # Only for one gene
gt(t(expr)[-c(50, 51), 1]~Exact_location * IBD*ID + AGE_SAMPLE + diagTime, data = meta_r[-c(50, 51), 1])

library("vegan") # For all lthe matrice
a <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime , meta_r[-c(50, 51), ], method = "euclidean")
# It captures the ~90% of the variation of the expression!!
a2 <- adonis2(as.data.frame(t(expr))[-c(50, 51), ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
a3 <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ Exact_location:ID , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
a4 <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ Exact_location:IBD , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
a5 <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ Exact_location , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
a6 <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ IBD , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
a7 <- adonis(as.data.frame(t(expr))[-c(50, 51), ] ~ AGE_SAMPLE , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")

b2 <- adonis(as.data.frame(t(otus_table_i))[-c(50, 51), ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
b3 <- adonis(as.data.frame(t(otus_table_i))[-c(50, 51), ] ~ Exact_location:ID , meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
b4 <- adonis(as.data.frame(t(otus_table_i))[-c(50, 51), ] ~ Exact_location:IBD, meta_r[-c(50, 51), ], method = "euclidean", by = "margin")
(b5 <- adonis(as.data.frame(t(otus_table_i))[-c(50, 51), ] ~ IBD , meta_r[-c(50, 51), ], method = "euclidean", by = "margin"))
