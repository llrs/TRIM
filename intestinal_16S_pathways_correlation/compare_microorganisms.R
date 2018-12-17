library("BioCor")
library("reactome.db")
library("org.Hs.eg.db")

today <- format(Sys.time(), "%Y%m%d")

#' Read correlation and creates data for BioCor
#' 
#' @param file Name of the file
#' @return A list of microorganisms to genes
read_cor <- function(file) {
  out <- read.csv(file)
  
  out <- split(out$Gene, out$Microorganism, drop = TRUE)
  lapply(out, function(x) {
    y <- as.character(x)
    unique(y[!is.na(y) & y != ""])
  })
  
}
cors_all <- readRDS("correlations_all.RDS")
cors_all <- read_cor("correlation_all.csv")
reactome <- as.list(reactomeEXTID2PATHID)
names(reactome) <- mapIds(org.Hs.eg.db, keys= names(reactome), keytype = "ENTREZID", column = "SYMBOL")
reactome <- reactome[!is.na(names(reactome))]

# Similarity between genes related to microorganisms
groups <- cors_all[lengths(cors_all) >= 10]
comp <- mclusterGeneSim(groups, reactome, c("BMA", "BMA"))
comp <- comp[is.finite(comp[1, ]), is.finite(comp[, 1])]

# MDS and plotting
mds <- cmdscale(as.dist(1-comp))
# 
png(paste0("Figures/", today, "_microorganisms_MDS.png"))
plot(mds*1.2, type = "n", xlab = "MDS1", ylab = "MDS2", main = "Microorganisms relation by function", asp = 1)
abline( h = 0, v = 0, col = "lightgrey")
text(mds[, 1], mds[, 2], labels = rownames(mds))
dev.off()
png(paste0("Figures/", today, "_microorganisms_MDS_zoom.png"))
sub_mds <- mds[abs(mds[, 1]) < 0.25 & abs(mds[, 2]) < 0.25, ]
plot(sub_mds, type = "n", xlab = "MDS1", ylab = "MDS2", main = "Microorganisms relation by function", asp = 1)
abline( h = 0, v = 0, col = "lightgrey")
text(sub_mds[, 1], sub_mds[, 2], labels = rownames(sub_mds))
dev.off()
