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

cors_all <- read_cor("correlation_all.csv")
reactome <- as.list(reactomeEXTID2PATHID)
names(reactome) <- mapIds(org.Hs.eg.db, keys= names(reactome), keytype = "ENTREZID", column = "SYMBOL")
reactome <- reactome[!is.na(names(reactome))]

# Similarity between genes related to microorganisms
groups <- cors_all[lengths(cors_all) >= 1]
comp <- mclusterGeneSim(groups, reactome, c("BMA", "BMA"))
comp <- comp[is.finite(comp[1, ]), is.finite(comp[, 1])]

# MDS and plotting
png(paste0("Figures/", today, "_microorganisms_MDS.png"))
mds <- cmdscale(comp)
plot(mds*1.1, type = "n", xlab = "MDS1", ylab = "MDS2", main = "Microorganisms relation by function")
text(mds[, 1], mds[, 2], labels = rownames(mds))
dev.off()
