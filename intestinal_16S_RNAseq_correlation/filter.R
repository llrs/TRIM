cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
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
expr <- as.matrix(read.delim(file.path(rna, "taula_sencera2.tsv"), 
                             check.names = FALSE))
file_meta_r <- file.path(rna, "metadata_28032018.csv")
meta_r <- read.table(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, sep = ";",
  na.strings = c(NA, ""), header = TRUE, dec = c(",", ".")
)

setwd(cd)

# Summarize to genus
library("metagenomeSeq")

# Create the objects to summarize data
MR_i <- newMRexperiment(
  otus_table_i,
  # phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")

# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- colnames(expr)[rev(position)]

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

# Subset if all the rows are 0 and if sd is 0
genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ]
expr <- expr[apply(expr, 1, sd) != 0, ]

abundance <- 0.005 # 0.5%

library("ComplexHeatmap")
library("org.Hs.eg.db")
library("heatmaply")

# Load data of correlations
all_s <- readRDS("correlations_all.RDS")
disease <- readRDS("correlations_CD.RDS")
controls <- readRDS("correlations_Controls.RDS")
disease_T0 <- readRDS("correlations_CD_T0.RDS")
disease_T26 <- readRDS("correlations_CD_T26.RDS")
disease_T52 <- readRDS("correlations_CD_T52.RDS")

# Load data of pvalues
pall_s <- readRDS("padj_all.RDS")
pdisease <- readRDS("padj_CD.RDS")
pcontrols <- readRDS("padj_Controls.RDS")
pdisease_T0 <- readRDS("padj_CD_T0.RDS")
pdisease_T26 <- readRDS("padj_CD_T26.RDS")
pdisease_T52 <- readRDS("padj_CD_T52.RDS")

filter_values <- function(file, cors, pval, threshold) {

  # Load data from all the patients
  pre <- "../intestinal_16S_RNAseq_metadb"
  load(file.path(pre, file))
  
  # Find outliers/important genes
  comp1 <- sgcca.centroid$a$RNAseq[, 1]
  # (RNAseq_sd <- sd(comp1))
  # (RNAseq_mean <- mean(comp1))
  # (RNAseq_median <- median(comp1))
  outliers <- comp1 != 0
  # summary(outliers)
  comp1 <- comp1[outliers]
  
  keepGenes <- colnames(cors) %in% names(comp1)
  cors <- cors[, keepGenes]
  pval <- pval[, keepGenes]
  
  colnames(cors) <- gsub("(.*)\\..*", "\\1", colnames(cors))
  
  symbol <- mapIds(
    org.Hs.eg.db, keys =  colnames(cors), keytype = "ENSEMBL",
    column = "SYMBOL"
  )
  colnames(cors) <- symbol
  colnames(pval) <- symbol
  cors <- cors[, !is.na(colnames(cors))]
  pval <- pval[, !is.na(colnames(pval))]
  
  message("Dimensions ", paste0(dim(cors), collapse = ", "))
  
  keepCols <- apply(pval, 2, function(x){any(x < threshold)})
  keepRows <- apply(pval, 1, function(x){any(x < threshold)})
  
  keepCols[is.na(keepCols)] <- FALSE
  keepRows[is.na(keepRows)] <- FALSE
  
  if (sum(keepCols) == 0 || sum(keepRows) == 0) {
    stop("No relevant correlations with this threshold")
  }
  
  cors <- cors[keepRows, keepCols]
  pval <- pval[keepRows, keepCols]
  if (is.null(pval) || is.null(cors)) {
    stop("No relevant correlations with this threshold")
  }
  message("Dimensions ", paste0(dim(cors), collapse = ", "))
  
  list(cors = cors, pval = pval)
}

relevant <- function(file, cors, pval, threshold = 0.05) {
  l <- filter_values(file, cors, pval, threshold)
  pval <- l$pval
  cors <- l$cors
  if (ncol(pval) == 0) {
    stop("No relevant correlations with this threshold")
  }
  ind <- as.data.frame(which(pval < threshold, arr.ind = TRUE), 
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind))
  cor_pval <- apply(ind, 1, function(x){
    c("cors" = cors[x[1], x[2]],
      "pvalue" = pval[x[1], x[2]])
  })
  ind$row <- rownames(cors)[ind$row]
  ind$col <- colnames(cors)[ind$col]
  ind <- cbind(ind, t(cor_pval))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")
  
  ind <- ind[!duplicated(ind), ]
  ind <- ind[order(ind$Microorganism, ind$pvalue, decreasing = c(TRUE, FALSE)), ]
  rownames(ind) <- seq_len(nrow(ind))
  ind
}

# Expects genes in rows and species at the columns
plot_cor <- function(file, cors, pval, threshold, label) {
  l <- filter_values(file, cors, pval, threshold)
  cors <- l$cors
  
  cors <- cors[!duplicated(rownames(cors)), ]
  colors_g <- ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                            midpoint = 0, limits = c(-1, 1))
  heatmaply(cors, name = "Cor",
            ylab = "Genes",
            xlab = "Genus",
            scale_fill_gradient_fun = colors_g,
            file = paste0("Figures/", today, "heatmap", label,".html"))

}

#  The numbers come from the number of samples in each correlation
threshold <- 0.05
b <- relevant("sgcca.RData", all_s, pall_s, threshold)
write.csv(b, file = "correlation_all.csv", row.names = FALSE, na = "")
# 
# b <- relevant("IBD.RData", disease, pdisease, threshold)
# write.csv(b, file = "correlation_CD.csv", row.names = FALSE, na = "")
# 
# b <- relevant("Controls.RData", controls, pcontrols, threshold)
# write.csv(b, file = "correlation_Controls.csv", row.names = FALSE, na = "")
# 
# b <- relevant("IBD_T0.RData", disease_T0, pdisease_T0, threshold)
# write.csv(b, file = "correlation_CD_T0.csv", row.names = FALSE, na = "")
# 
# b <- relevant("IBD_T26.RData", disease_T26, pdisease_T26, threshold)
# write.csv(b, file = "correlation_CD_T26.csv", row.names = FALSE, na = "")
# 
# b <- relevant("IBD_T52.RData", disease_T52, pdisease_T52, threshold)
# write.csv(b, file = "correlation_CD_T52.csv", row.names = FALSE, na = "")

plot_single_cor <- function(x, gene, y, rowY, colr, case) {
  suppressMessages(g <- mapIds(org.Hs.eg.db, keys = gene, 
                               keytype = "SYMBOL", column = "ENSEMBL"))
  rowX <- grep(g, rownames(x))
  x_s <- x[rowX, ]
  x_s[x_s == 0] <- NA
  
  y_s <- y[rowY, ]
  y_s[y_s == 0] <- NA
  
  x_s <- log10(x_s)
  y_s <- log10(y_s)
  main <- cor(x_s, y_s, method = "spearman", use = "pairwise.complete.obs")
  plot(x_s, y_s, xlab = gene, ylab = rowY, main = main, col = colr, 
       pch = as.numeric(case))
  legend("topright", fill = as.factor(levels(colr)), legend = levels(colr))
  legend("bottomleft", pch = as.factor(levels(case)), 
         legend = levels(case))
  
}

pdf("correlations_plot.pdf")
apply(b, 1, function(x) {
  micro <- x[[1]]
  gene <- x[[2]]
  plot_single_cor(expr, gene, genus_i, micro, as.factor(meta_r$IBD), as.factor(meta_r$IBD))
})
dev.off()