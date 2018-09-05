# Load ####
cd <- setwd("..")

# Load genes of chromosome X and Y
source("../genes_XY.R")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

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
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(cd)

# Summarize to genus ####
library("metagenomeSeq")
library("vegan")
library("phyloseq")

# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# Create the objects to summarize data
MR_i <- newMRexperiment(
  otus_table_i,
  # phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix", norm = TRUE, log = TRUE)

colnames(genus_i) <- gsub("[0-9]+\\.(.*)", "\\1", colnames(genus_i))

gsv <- readRDS("../intestinal_16S_pathways_metadb/GSV.RDS")
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(genus_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(gsv), ]

# Subset expression and outs
genus_i <- genus_i[, meta_r$Seq_code_uDNA]
gsv <- gsv[, meta_r$`Sample Name_RNA`]

# Filter by those which have more than one data point
genus_i <- genus_i[apply(genus_i, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
gsv <- gsv[apply(gsv, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]

# Filter by abundance at 0.5%
abundance <- 0.005 # 0.5%
a <- prop.table(genus_i, 2)
b <- rowSums(a > abundance)
genus_i <- genus_i[b != 0, ]

saveRDS(gsv, "gsv.RDS")
saveRDS(genus_i, "genus.RDS")

## Functions ####
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
  
  # Filter by sex
  sexual_related <- gsub("(.+)\\..*", "\\1", rownames(expr)) %in% 
    c(#bmX$ensembl_gene_id, 
      bmY$ensembl_gene_id)
  expr <- expr[!sexual_related, ]
  
  expr[expr == 0] <- NA
  genus_i[genus_i == 0] <- NA

  r <- matrix(NA, nrow = nrow(genus_i), ncol = nrow(expr), 
              dimnames = list(rownames(genus_i), rownames(expr)))
  pval <- matrix(NA, nrow = nrow(genus_i), ncol = nrow(expr), 
                 dimnames = list(rownames(genus_i), rownames(expr)))
  
  for (gene in rownames(expr)) {
    if (sum(!is.na(y)) < 3 | sum(!is.na(y))/length(y) < 0.15) {
      next
    }
    gene_expr <- expr[gene, ]
    for (micro in rownames(genus_i)) {
      if (sum(!is.na(x)) < 3 | sum(!is.na(x))/length(x) < 0.15) {
        next
      }
      genus_expr <- genus_i[micro, ]
      # Join all the data
      d <- rbind(gene_expr, genus_expr)
      pairwise <- apply(d, 2, function(z){all(!is.na(z))})
      k <- sum(pairwise, na.rm = TRUE)
      if (k/length(pairwise) < 0.15 & k < 4) {
        next
      }
      
      # Errors because they don't match up to three points
      try({
        cors <- cor.test(gene_expr, genus_expr, use = "spearman", 
                         use = "pairwise.complete.obs")
        pval[micro, gene] <- cors$p.value
        r[micro, gene] <- cors$estimate},
        silent = TRUE)
    }
  }
  # pval <- pvalue(p, ncol(expr))
  padj <- p.adjust(pval, method = "BH")
  dim(padj) <- dim(pval)
  dimnames(padj) <- dimnames(pval)
  saveRDS(r, file = paste0("correlations_", label, ".RDS"))
  saveRDS(padj, file = paste0("padj_", label, ".RDS"))
}

## All samples ####
# ncol(expr)
cor2(genus_i, gsv, "all")


## CD ####
keep <- meta_r$IBD == "CD"
# sum(keep)
# cor_sign(sum(keep))
cor2(genus_i[, keep], gsv[, keep], "CD")

## Healthy ####
keep <-  meta_r$IBD == "CONTROL"
# sum(keep)
# cor_sign(sum(keep))
cor2(genus_i[, keep], gsv[, keep], "Controls")
