cd <- setwd("..")

source("genes_XY.R") # for BMY object
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


# Normalize expression
norm_edgeR <- function(expr){
  expr_edge <- edgeR::DGEList(expr)
  expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
  edgeR::cpm(expr_edge, normalized.lib.sizes = TRUE, log = TRUE)
}

# Filter expression
expr_norm <- norm_RNAseq(expr)
expr <- norm_RNAseq(expr_norm)

# Create the objects to summarize data
MR_i <- newMRexperiment(
  otus_table_i,
  # phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)

MR_in <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
species_i <- aggTax(MR_in, lvl = "Species", out = "matrix", norm = TRUE, log = TRUE)
genus_i <- aggTax(MR_in, lvl = "Genus", out = "matrix", norm = TRUE, log = TRUE)

# Subset if all the rows are 0 and if sd is 0
species_i <- species_i[apply(species_i, 1, sd) != 0, ]
genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ]
expr <- expr[apply(expr, 1, sd) != 0, ]

# Filter by those which have more than one data point
species_i <- species_i[apply(species_i, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
genus_i <- genus_i[apply(genus_i, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
expr <- expr[apply(expr, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]

# Filter by abundance at 0.5%
species_ie <- aggTax(MR_i, lvl = "Species", out = "matrix", norm = FALSE, log = FALSE)
genus_ie <- aggTax(MR_i, lvl = "Genus", out = "matrix", norm = FALSE, log = FALSE)

# Subset if all the rows are 0 and if sd is 0
species_ie <- species_ie[apply(species_ie, 1, sd) != 0, ]
genus_ie <- genus_ie[apply(genus_ie, 1, sd) != 0, ]

# Filter by those which have more than one data point
species_ie <- species_ie[apply(species_ie, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]
genus_ie <- genus_ie[apply(genus_ie, 1, function(x) {sum(x != 0)/length(x) > 0.15}), ]

abundance <- 0.005 # 0.5%
a <- prop.table(genus_ie, 2)
b <- rowSums(a > abundance)
genus_i <- genus_i[b != 0, ]

a <- prop.table(species_ie, 2)
b <- rowSums(a > abundance)
species_i <- species_i[b != 0, ]

saveRDS(expr, "expr.RDS")
saveRDS(species_i, "species.RDS")
saveRDS(genus_i, "genus.RDS")

## Functions ####
#' Correlation matrix
#' 
#' Calculates the correlation between species at 0.005
#' @param x The species table
#' @param y the expression matrix
#' @param label The label of the correlation
#' @param abundance The proportion below to filter out
#' @return A file of correlations
cor2 <- function(x, y, label, abundance = 0.005){
  species_i <- as.matrix(x)
  expr <- as.matrix(y)
  
  # Filter by sex
  sexual_related <- trimVer(rownames(expr)) %in% 
    c(#bmX$ensembl_gene_id, 
      bmY$ensembl_gene_id)
  expr <- expr[!sexual_related, ]
  
  expr[expr == 0] <- NA
  species_i[species_i == 0] <- NA

  r <- matrix(NA, nrow = nrow(species_i), ncol = nrow(expr), 
              dimnames = list(rownames(species_i), rownames(expr)))
  pval <- matrix(NA, nrow = nrow(species_i), ncol = nrow(expr), 
                 dimnames = list(rownames(species_i), rownames(expr)))
  
  for (gene in rownames(expr)) {
    gene_expr <- expr[gene, ]
    if (sum(!is.na(gene_expr)) < 3 | sum(!is.na(gene_expr))/length(gene_expr) < 0.15) {
      next
    }
    for (micro in rownames(species_i)) {
      species_expr <- species_i[micro, ]
      if (sum(!is.na(species_expr)) < 3 | sum(!is.na(species_expr))/length(species_expr) < 0.15) {
        next
      }
      # Join all the data
      d <- rbind(gene_expr, species_expr)
      pairwise <- apply(d, 2, function(z){all(!is.na(z))})
      k <- sum(pairwise, na.rm = TRUE)
      if (k/length(pairwise) < 0.15 & k < 4) {
        next
      }
      
      # Errors because they don't match up to three points
      try({
        cors <- cor.test(gene_expr, species_expr, use = "spearman", 
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
  saveRDS(pval, file = paste0("pval_", label, ".RDS"))
}

## All samples ####
# ncol(expr)
cor2(species_i, expr, "all_species")
cor2(genus_i, expr, "all_genus")

## CD ####
keep <- meta_r$IBD == "CD"
cor2(genus_i[, keep], expr[, keep], "CD_genus")
cor2(species_i[, keep], expr[, keep], "CD_species")

## Healthy ####
keep <-  meta_r$IBD == "CONTROL"
cor2(genus_i[, keep], expr[, keep], "Controls_genus")
cor2(species_i[, keep], expr[, keep], "Controls_species")

# Location ####
ileum <-  meta_r$Exact_location %in% "ILEUM"
cor2(genus_i[, ileum], expr[, ileum], "Ileum_genus")
cor2(species_i[, ileum], expr[, ileum], "Ileum_species")

cor2(genus_i[, !ileum], expr[, !ileum], "Colon_genus")
cor2(species_i[, !ileum], expr[, !ileum], "Colon_species")
