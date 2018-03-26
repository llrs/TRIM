cd <- setwd("..")

# Load the helper file
library("integration")
library("ppcor")

intestinal <- "intestinal_16S"

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)

tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Load the input data
rna <- "intestinal_RNAseq"
expr <- as.matrix(read.delim(file.path(rna, "table.counts.results"), 
                             check.names = FALSE))
file_meta_r <- file.path(rna, "metadata_13032018.csv")
meta_r <- read.table(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, sep = ";",
  na.strings = c(NA, ""), header = TRUE, dec = c(",", ".")
)

setwd(cd)

# Summarize to genus
library("metagenomeSeq")
library("vegan")
library("phyloseq")

# Create the objects to summarize data
MR_i <- newMRexperiment(
  otus_table_i,
  # phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")

# Calculate the alpha diversity of the samples
alpha <- estimate_richness(otu_table(genus_i, taxa_are_rows = TRUE))
alpha$Shannon.Effective <- exp(alpha$Shannon)
# Remove some diversities that are not linear


# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- rev(position)

# Find the samples that we have microbiota and expression
int <- intersect(
  meta_r$Sample_Code_uDNA[!is.na(meta_r$Sample_Code_uDNA) &
    !is.na(meta_r$`Sample Name_RNA`)],
  meta_i$Sample_Code
)

meta_i <- meta_i[meta_i$Sample_Code %in% int, ]
meta_r <- meta_r[meta_r$Sample_Code_uDNA %in% int, ]
meta_r <- meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Match the labels and order to append the id
meta_i <- meta_i[match(meta_r$Sample_Code_uDNA, meta_i$Sample_Code), ]
meta_r$`Sample Name_Code` <- gsub("([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1", rownames(meta_i))

colnames(genus_i) <- gsub(
  "([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1",
  colnames(genus_i)
)

# Subset expression and outs
expr <- expr[, meta_r$`Sample Name_RNA`]
genus_i <- genus_i[, meta_r$`Sample Name_Code`]
alpha <- alpha[paste0("X", meta_r$`Sample Name_Code`), ]

# Subset if all the rows are 0 and if sd is 0
genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ]
expr <- expr[apply(expr, 1, sd) != 0, ]

abundance <- 0.005 # 0.5%

cors <- cor(log10(t(expr) + 0.25), alpha)

## All samples ####
#' Correlation matrix
#' 
#' Calculates the correlation between genus at 0.005
#' @param x The genus table
#' @param y the expression matrix
#' @param label The label of the correlation
#' @param abundance The proportion below to filter out
#' @return A file of correlations
cor2 <- function(x, y, label, abundance = 0.005){
  genus_i <- as.matrix(x)
  expr <- as.matrix(y)
  
  # Filter by those which have variation
  genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ]
  expr <- expr[apply(expr, 1, sd) != 0, ]
  
  # Filter by those which have more than one data point
  genus_i <- genus_i[apply(genus_i, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
  expr <- expr[apply(expr, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
  
    # Filter by abundance at 0.5%
  a <- prop.table(genus_i, 2)
  b <- rowSums(a > abundance)
  
  genus_i <- genus_i[b != 0, ]
  
  # Correlate
  p <- cor(log10(t(expr +0.25)), log10(t(genus_i+0.25)), method = "spearman")
  # Subset to only between microbiota and RNAseq
  # p$estimates <- p$estimates[rownames(genus_i), rownames(expr)]
  # p$p.value <- p$p.value[rownames(genus_i), rownames(expr)]
  # p$statistic <- p$statistic[rownames(genus_i), rownames(expr)]
  pval <- pvalue(p, ncol(expr))
  padj <- p.adjust(pval, method = "BH")
  dim(padj) <- dim(pval)
  dimnames(padj) <- dimnames(pval)
  saveRDS(p, file = paste0("correlations_", label, ".RDS"))
  saveRDS(padj, file = paste0("padj_", label, ".RDS"))
}

ncol(expr)
cor2(genus_i, expr, "all")

keep <- meta_r$IBD == "CD"
sum(keep)
cor_sign(sum(keep))
cor2(genus_i[, keep], expr[, keep], "CD")

keep <-  meta_r$IBD == "CONTROL"
sum(keep)
cor_sign(sum(keep))
cor2(genus_i[, keep], expr[, keep], "Controls")

keep <-  meta_r$IBD == "CD" & meta_r$Time == "T0"
sum(keep)
cor_sign(sum(keep))
cor2(genus_i[, keep], expr[, keep], "CD_T0")

keep <-  meta_r$IBD == "CD" & meta_r$Time == "T26"
sum(keep)
cor_sign(sum(keep))
cor2(genus_i[, keep], expr[, keep], "CD_T26")

keep <-  meta_r$IBD == "CD" & meta_r$Time == "T52"
sum(keep)
cor_sign(sum(keep))
cor2(genus_i[, keep], expr[, keep], "CD_T52")
