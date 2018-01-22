
source("../helper_functions.R")

# Load data of correlations
cors <- readRDS("correlations.RDS")
# Load data from all the patients
pre <- "../intestinal_16S_RNAseq_metadb"
load(file.path(pre, "sgcca.RData"))

file_meta_r <- "../intestinal_RNAseq/111217_metadata.csv"
meta_r <- read.table(file_meta_r, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = ";",
                     na.strings = c(NA, ""), header = TRUE, dec = c(",", "."))


# Correct metadata
meta_r <- meta_r_norm(meta_r)

# Find outliers
comp1 <- sgcca.centroid$a$RNAseq[, 1]
(RNAseq_sd <- sd(comp1))
(RNAseq_mean <- mean(comp1))
(RNAseq_median <- median(comp1))
outliers <- comp1 != 0

library("org.Hs.eg.db")
names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))

cors <- cors[, names(comp1)]
symbol <- mapIds(org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL", 
                   column = "SYMBOL")

outliers[is.na(symbol)] <- FALSE

# pdf(paste0("Figures/", today, "_plots.pdf"))
gplots::heatmap.2(cors[, outliers], main = "Correlation heatmap: genes-genus", 
                  xlab = "Genes", ylab = "Genus", scale = "none", 
                  tracecol = "black", col = bluered(11), trace = "none", 
                  labCol = symbol[outliers], 
                  colCol = ifelse(comp1[outliers] > 0, "green","black"))
