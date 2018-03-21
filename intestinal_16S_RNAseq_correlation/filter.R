
library("integration")
library("ComplexHeatmap")
library("org.Hs.eg.db")

# Load data of correlations
all_s <- readRDS("correlations_all.RDS")
disease <- readRDS("correlations_CD.RDS")
controls <- readRDS("correlations_Controls.RDS")
disease_T0 <- readRDS("correlations_CD_T0.RDS")
disease_T26 <- readRDS("correlations_CD_T26.RDS")
disease_T52 <- readRDS("correlations_CD_T52.RDS")



plot_cor <- function(file, cors, n) {
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
  cors <- cors[, names(comp1)]
  
  names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))
  
  symbol <- mapIds(
    org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL",
    column = "SYMBOL"
  )
  colnames(cors) <- symbol
  mat <- cors[, !is.na(colnames(cors))]
  
  message("Dimensions ", paste0(dim(cors), collapse = ", "))
  
  d2 <- apply(mat, 1:2, function(x){
    p.adjust(pvalue(x, n = n), n = ncol(mat)*nrow(mat))
  })
  keepCols <- apply(d2, 2, function(x){any(x < 0.05)})
  keepRows <- apply(d2, 1, function(x){any(x < 0.05)})
  d2 <- d2[keepRows, keepCols]
  mat <- mat[keepRows, keepCols]
  
  message("Dimensions ", paste0(dim(mat), collapse = ", "))
  print(Heatmap(mat, name = "Cor", 
                column_title = "Genes",
                row_title = "Genus"))
  
  ind <- as.data.frame(which(d2 < 0.05, arr.ind = TRUE), 
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind))
  # stopifnot(length(ind$row) != length(rownames(mat)))
  # stopifnot(length(ind$row) != length(colnames(mat)))
  cors <- apply(ind, 1, function(x){
    c("cors" = mat[x[1], x[2]],
      "pvalue" = d2[x[1], x[2]])
  })
  ind$row <- rownames(mat)[ind$row]
  ind$col <- colnames(mat)[ind$col]
  ind <- cbind(ind, t(cors))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")
  
  ind <- ind[!duplicated(ind), ]
  rownames(ind) <- seq_len(nrow(ind))
  return(invisible(ind))
}

#  The numbers come from the number of samples in each correlation
#  
pdf(paste0("Figures/", today, "_plots.pdf"))
b <- plot_cor("sgcca.RData", all_s, 152)
write.csv(b, file = "correlation_all.csv", row.names = FALSE, na = "")

b <- plot_cor("IBD.RData", disease, 101)
write.csv(b, file = "correlation_CD.csv", row.names = FALSE, na = "")

b <- plot_cor("Controls.RData", controls, 51)
write.csv(b, file = "correlation_Controls.csv", row.names = FALSE, na = "")

b <- plot_cor("IBD_T0.RData", disease_T0, 15)
write.csv(b, file = "correlation_CD_T0.csv", row.names = FALSE, na = "")

b <- plot_cor("IBD_T26.RData", disease_T26, 28)
write.csv(b, file = "correlation_CD_T26.csv", row.names = FALSE, na = "")

b <- plot_cor("IBD_T52.RData", disease_T52, 26)
write.csv(b, file = "correlation_CD_T52.csv", row.names = FALSE, na = "")

dev.off()

# To see if the weight have a relation with the direct correlation
# colCol = ifelse(comp1[outliers] > 0, "green", "black")

# To have a good feeling plot some random genes. Most of the time there are
# the double of 0 than in the selected heatmap
# sam <- sample(colnames(cors), 600)
# gplots::heatmap.2(cors[, sam], main = "Correlation heatmap: genes-genus",
#                   xlab = "Genes", ylab = "Genus", scale = "none",
#                   tracecol = "black", col = bluered(64), trace = "none",
#                   labCol = symbol[sam], margins = c(6, 9))
# 
# disease <- disease[, colnames(disease) %in% names(outliers)]
# a <- matrix(, ncol = ncol(disease[, outliers]), nrow = nrow(cors))
# a[abs(disease[, outliers]) >= 0.19] <- "*" # Significant threshold of 0.05
# heatmap.2(
#   disease[, outliers], main = "Correlation heatmap IBD: genes-genus",
#   xlab = "Genes", ylab = "Genus", scale = "none",
#   tracecol = "black", col = bluered(64), trace = "none",
#   labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
#   notecol = "black", notecex = 0.5
# )
# 
# 
# controls <- controls[, colnames(controls) %in% names(outliers)]
# a <- matrix(, ncol = ncol(controls[, outliers]), nrow = nrow(cors))
# a[abs(controls[, outliers]) >= 0.28] <- "*" # Significant threshold of 0.05
# heatmap.2(
#   controls[, outliers], main = "Correlation heatmap controls: genes-genus",
#   xlab = "Genes", ylab = "Genus", scale = "none",
#   tracecol = "black", col = bluered(64), trace = "none",
#   labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
#   notecol = "black", notecex = 0.5
# )
# 
# # Do the same but using the genes of the model with only IBD
# pre <- "../intestinal_16S_RNAseq_metadb"
# load(file.path(pre, "IBD.RData"))
# comp1 <- sgcca.centroid$a[["RNAseq"]][, 1]
# outliers <- comp1 != 0
# summary(outliers)
# 
# cors <- readRDS("correlations.RDS")
# disease <- readRDS("correlations_IBD.RDS")
# controls <- readRDS("correlations_C.RDS")
# 
# cors <- cors[, names(comp1)]
# 
# library("org.Hs.eg.db")
# names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))
# 
# symbol <- mapIds(
#   org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL",
#   column = "SYMBOL"
# )
# 
# outliers[is.na(symbol)] <- FALSE
# 
# 
# a <- matrix(, ncol = ncol(cors[, outliers]), nrow = nrow(cors))
# a[abs(cors[, outliers]) >= 0.16] <- "*" # Significant threshold of 0.05
# heatmap.2(
#   cors[, outliers], main = "Correlation heatmap all: genes-genus",
#   xlab = "Genes", ylab = "Genus", scale = "none",
#   tracecol = "black", col = bluered(64), trace = "none",
#   labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
#   notecol = "black", notecex = 0.5
# )
# 
# disease <- disease[, colnames(disease) %in% names(outliers)]
# a <- matrix(, ncol = ncol(disease[, outliers]), nrow = nrow(cors))
# a[abs(disease[, outliers]) >= 0.19] <- "*" # Significant threshold of 0.05
# heatmap.2(
#   disease[, outliers], main = "Correlation heatmap IBD: genes-genus",
#   xlab = "Genes", ylab = "Genus", scale = "none",
#   tracecol = "black", col = bluered(64), trace = "none",
#   labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
#   notecol = "black", notecex = 0.5
# )
# 
# int <- intersect(colnames(controls), names(outliers))
# controls <- controls[, int]
# a <- matrix(, ncol = ncol(controls), nrow = nrow(cors))
# a[abs(controls[, outliers[int]]) >= 0.28] <- "*" # Significant threshold of 0.05
# heatmap.2(
#   controls[, outliers[int]], main = "Correlation heatmap controls: genes-genus",
#   xlab = "Genes", ylab = "Genus", scale = "none",
#   tracecol = "black", col = bluered(64), trace = "none",
#   labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
#   notecol = "black", notecex = 0.5
# )
