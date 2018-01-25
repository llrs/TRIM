
source("../helper_functions.R")

# Load data of correlations
cors <- readRDS("correlations.RDS")
disease <- readRDS("correlations_IBD.RDS")
controls <- readRDS("correlations_C.RDS")

# Load data from all the patients
pre <- "../intestinal_16S_RNAseq_metadb"
load(file.path(pre, "sgcca.RData"))

# Find outliers/important genes
comp1 <- sgcca.centroid$a$RNAseq[, 1]
(RNAseq_sd <- sd(comp1))
(RNAseq_mean <- mean(comp1))
(RNAseq_median <- median(comp1))
outliers <- comp1 != 0
summary(outliers)

cors <- cors[, names(comp1)]

library("org.Hs.eg.db")
names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))

symbol <- mapIds(
  org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL",
  column = "SYMBOL"
)

outliers[is.na(symbol)] <- FALSE

pdf(paste0("Figures/", today, "_plots.pdf"))
library("gplots")
a <- matrix(, ncol = ncol(cors[, outliers]), nrow = nrow(cors))
a[abs(cors[, outliers]) >= 0.16] <- "*" # Significant threshold of 0.05
heatmap.2(
  cors[, outliers], main = "Correlation heatmap all: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)
# To see if the weight have a relation with the direct correlation
# colCol = ifelse(comp1[outliers] > 0, "green", "black")

# To have a good feeling plot some random genes. Most of the time there are
# the double of 0 than in the selected heatmap
# sam <- sample(colnames(cors), 600)
# gplots::heatmap.2(cors[, sam], main = "Correlation heatmap: genes-genus",
#                   xlab = "Genes", ylab = "Genus", scale = "none",
#                   tracecol = "black", col = bluered(64), trace = "none",
#                   labCol = symbol[sam], margins = c(6, 9))

disease <- disease[, colnames(disease) %in% names(outliers)]
a <- matrix(, ncol = ncol(disease[, outliers]), nrow = nrow(cors))
a[abs(disease[, outliers]) >= 0.19] <- "*" # Significant threshold of 0.05
heatmap.2(
  disease[, outliers], main = "Correlation heatmap IBD: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)


controls <- controls[, colnames(controls) %in% names(outliers)]
a <- matrix(, ncol = ncol(controls[, outliers]), nrow = nrow(cors))
a[abs(controls[, outliers]) >= 0.28] <- "*" # Significant threshold of 0.05
heatmap.2(
  controls[, outliers], main = "Correlation heatmap controls: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)
