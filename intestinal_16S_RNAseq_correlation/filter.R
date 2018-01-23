
source("../helper_functions.R")

# Load data of correlations
cors <- readRDS("correlations.RDS")
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

symbol <- mapIds(org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL", 
                   column = "SYMBOL")

outliers[is.na(symbol)] <- FALSE

pdf(paste0("Figures/", today, "_plots.pdf"))
library("gplots")
heatmap.2(cors[, outliers], main = "Correlation heatmap: genes-genus", 
                  xlab = "Genes", ylab = "Genus", scale = "none", 
                  tracecol = "black", col = bluered(64), trace = "none", 
                  labCol = symbol[outliers], margins = c(6, 9))
    # To see if the weight have a relation with the direct correlation
                  # colCol = ifelse(comp1[outliers] > 0, "green", "black") 

# To have a good feeling plot some random genes. Most of the time there are 
# the double of 0 than in the selected heatmap
# sam <- sample(colnames(cors), 600)
# gplots::heatmap.2(cors[, sam], main = "Correlation heatmap: genes-genus", 
#                   xlab = "Genes", ylab = "Genus", scale = "none", 
#                   tracecol = "black", col = bluered(64), trace = "none", 
#                   labCol = symbol[sam], margins = c(6, 9))
