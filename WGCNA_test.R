library("WGCNA")
enableWGCNAThreads()
data.wgcna <- t(otus_table_i[meta_i$HSCT_responder != "C", ])
library("phyloseq")
data("GlobalPatterns")
data.wgcna <- t(otu_table(GlobalPatterns))
powers <- seq(from = 1, to = 20, by = 1)
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, networkType = "unsigned", RsquaredCut = 0.8, corFnc = "bicor")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1], sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
  main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

net <- blockwiseModules(
  apply(data.wgcna, 1:2, as.numeric), power = sft$powerEstimate + 1,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  saveTOMFileBase = "stools_16S_TOM",
  verbose = 3, corType = "bicor"
)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)

table(mergedColors)
