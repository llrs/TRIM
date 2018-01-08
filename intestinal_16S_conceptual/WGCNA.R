wd <- setwd("..")
source("helper_functions.R")
intestinal <- "intestinal_16S"

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
setwd(wd)


# Clean the metadata
meta_i <- meta_i[, apply(meta_i, 2, function(x){length(unique(x)) != 1})]
meta_i$ID <- meta_i$Patient_ID
meta_i$ID[meta_i$Patient_ID %in% c("15", "23")] <- "15/23"
meta_i$ID[meta_i$Patient_ID %in% c("33", "36")] <- "33/36"
meta_i$ID[meta_i$Patient_ID %in% c("29", "35")] <- "29/35"
meta_i$ID <- as.factor(meta_i$ID)
# There is a mislabeling on those tubes, we don't know which is which
meta_i$CD_Aftected_area[meta_i$Sample_Code == "22_T52_T_DM_III"] <- NA


library("WGCNA")
enableWGCNAThreads()
# All together is not possible
data.wgcna <- t(otus_table_i[, meta_i$HSCT_responder != "C" & meta_i$CD_Aftected_area %in% "COLON"])


powers <-  seq(from = 1, to=20, by=1)
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, networkType = "unsigned", RsquaredCut = 0.8)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 8
adjacency = adjacency(data.wgcna, power = softPower, type = "unsigned");

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 5;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEs = moduleEigengenes(data.wgcna, colors = dynamicColors)$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
