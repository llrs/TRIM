library("integration")
library("org.Hs.eg.db")
library("fastDummies")
library("WGCNA")
library("RGCCA")
library("pheatmap")
library("gplots")
library("WriteXLS")

model_RGCCA <- readRDS("../intestinal_16S_RNAseq_metadb/model3_forced_interaction.RDS")

genes <- abs(sort(model_RGCCA$a[[1]][, 1]))
micros <- abs(sort(model_RGCCA$a[[2]][, 1]))
n_genes <- length(genes)
n_micros <- length(micros)
vars <- unlist(sapply(model_RGCCA$a[3:5], function(x){x[, 1]}))

non_empty <- function(x) {
    x[x != 0]
}
genes <- non_empty(genes)
micros <- non_empty(micros)

# Plot ecdf weights ####
plot(ecdf(micros), main = "ECDF", xlab = "weight", pch = NA)
plot(ecdf(genes), add = TRUE, col = "darkgreen")
legend("bottomright", fill = c("black", "darkgreen"), legend = c("OTUs", "Genes"))

# Look for number of correlations ####
q_threshold <- 0.95
rel_micros <- micros[micros > quantile(micros, q_threshold)]
rel_genes <- genes[genes > quantile(genes, q_threshold)]
q <- seq(from = 0, to = 0.95, by = 0.05)
pall_s <- readRDS("../intestinal_16S_RNAseq_correlation/padj_all_genus.RDS")
otus_tax <- readRDS("../intestinal_16S_RNAseq_metadb/otus_tax.RDS")
df <- sapply(q, function(x) {
    k_genes <- genes > quantile(genes, x)
    k_micros <- micros > quantile(micros, x)
    k_vars <- vars > quantile(vars, x)
    s_micros <- rownames(pall_s) %in% otus_tax[names(k_micros)[k_micros], "Genus"]
    s_genes <- colnames(pall_s) %in% names(genes)[k_genes]
    p001 <- sum(pall_s[s_micros, s_genes] < 0.001)
    p01 <- sum(pall_s[s_micros, s_genes] < 0.01)
    p05 <- sum(pall_s[s_micros, s_genes] < 0.05)
    p1 <- sum(pall_s[s_micros, s_genes] < 0.1)
    c(q = 1 - x,
      n_micros = sum(k_micros),
      n_genes = sum(k_genes),
      n_vars = sum(k_vars),
      sign_0.001 = p001,
      sign_0.01 = p01,
      sign_0.05 = p05,
      sign_0.1 = p1)
})
(df <- as.data.frame(t(df)))

s_micros <- rownames(pall_s) %in% otus_tax[names(micros), "Genus"]
s_genes <- colnames(pall_s) %in% names(genes)
p_threshold <- 0.001
pos <- which(pall_s[s_micros, s_genes] < p_threshold, arr.ind = TRUE)


org <- rownames(pos)
pos <- as.data.frame(pos)
pos$genes <- colnames(pall_s[s_micros, s_genes])[pos$col]
pos$org <- org

## Test if there are more correlations than expected ####
t1 <- table(pall_s[s_micros, s_genes] < p_threshold)
t2 <- table(pall_s[!s_micros, !s_genes] < p_threshold)

ts <- rbind(t1, t2)
fisher.test(ts) # p-value = 0.0008831 to 0.0008573
chisq.test(ts)  # p-value = 0.002665 to 0.002763

subsets_genes <- lapply(q, function(x) {
    k_genes <- genes > quantile(genes, x)
    names(genes[k_genes])
})

names(subsets_genes) <- paste("quantile", 1 - q, sep = "_")

# fgsea::fgsea(subset_genes, cor_genes)

k <- mapIds(org.Hs.eg.db, keys = trimVer(names(rel_genes)), keytype = "ENSEMBL", column = "SYMBOL")

meta_r <- readRDS("../intestinal_16S_RNAseq_metadb/meta.RDS")

# Plot 2nd dimension  ####

location <- ifelse(meta_r$Exact_location == "ILEUM", "Ileum", "Colon")
# Make the dimension plot
plot(model_RGCCA$Y[[1]][, 2], model_RGCCA$Y[[2]][, 2], xlab = "transcriptome", ylab = "microbiome",
    main = "Second dimensions", col = as.factor(meta_r$IBD), pch = 16)
legend("bottomleft", fill = c("black", "red"), legend = c("CD", "Controls"))
plot(model_RGCCA$Y[[1]][, 2], model_RGCCA$Y[[2]][, 2], xlab = "transcriptome", ylab = "microbiome",
    main = "Second dimensions", col = as.factor(location), pch = 16)
legend("bottomleft", fill = c("black", "red"), legend = c("Ileum", "Colon"))

# Recalculate the dummy variables for the names
m3 <- c("ID","SEX", "Surgery", "Treatment")
A3_2 <- dummy_cols(meta_r[, m3], m3, remove_first_dummy = TRUE)
A3cols <- colnames(A3_2[, !colnames(A3_2) %in% m3])

m5 <- "Transplant"
A5_2 <- dummy_cols(meta_r[, m5, drop = FALSE], m5, remove_first_dummy = TRUE)
A5cols <- c("AgeDiag", "AGE_SAMPLE", colnames(A5_2[, !colnames(A5_2) %in% m5]))

m4 <- c("Exact_location")
A4_2 <- dummy_cols(meta_r[, m4, drop = FALSE], m4, remove_first_dummy = TRUE)
A4cols <- colnames(A4_2[, !colnames(A4_2) %in% m4])

pdf("Figures/20190614_weights_model_RGCCA.pdf")
plot(c(-0.6, 1), c(-1, 0.6),  type = "n", ylab = "Dim 2", xlab = "Dim 1", main = "Weights")
keep1 <- rowSums(model_RGCCA$a[[1]]) != 0
keep2 <- rowSums(model_RGCCA$a[[2]]) != 0
points(model_RGCCA$a[[1]][keep1, ], pch = 15, col = "green")
points(model_RGCCA$a[[2]][keep2, ], pch = 17, col = "red")
# points(model_RGCCA$a[[3]], pch = 16, col = "brown")
text(model_RGCCA$a[[3]], labels = A3cols)
# points(model_RGCCA$a[[4]], pch = 16, col = "brown")
text(model_RGCCA$a[[4]], labels = A4cols)
# points(model_RGCCA$a[[5]], pch = 16, col = "brown")
text(model_RGCCA$a[[5]], labels = A5cols)
legend("topright", fill = c("green", "red", "brown"), legend = c("RNA", "DNA", "Clinical variables"))
dev.off()

# Heatmap ####

w <- do.call(rbind, model_RGCCA$a)
r <- c(rownames(model_RGCCA$a[[1]]), rownames(model_RGCCA$a[[2]]), A3cols, A4cols, A5cols)
rownames(w) <- r
empty <- apply(w, 1, function(x){all(x == 0)})
w <- w[!empty, ]

d <- dist(w)
dm <- as.matrix(d)
heatmap(dm[7500:7745, 7500:7745])

# WGCNA ####

## Data normalized separately ####
A <- readRDS("../intestinal_16S_RNAseq_metadb/model3_TRIM.RDS")
bigmatrix <- cbind(A[[1]][, names(genes)], A[[2]][, names(micros)])
biggermatrix <- cbind(A[[1]], A[[2]])

## Apply WGCNA ####
allowWGCNAThreads(3)
input <- bigmatrix
sdata <- apply(input, 2, scale2)
rownames(sdata) <- rownames(input)
fulldata <- t(sdata)
colnames(fulldata) <- rownames(sdata)

powers <- seq(1, 20, by = 1)
psoft <- pickSoftThreshold(data = sdata, powerVector = powers, verbose = TRUE)
sft <- psoft

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
par(mfrow = c(1, 1))

# taking power=7 remembering previous info prom LLuis
adjacency <- adjacency(sdata, power = 7)
dissTOM <- 1 - TOMsimilarity(adjacency)
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 10
dynamicMods <- cutreeDynamic(dendro = geneTree,  method = "tree", 
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
scores <- rbind(model_RGCCA$a[[1]], model_RGCCA$a[[2]])
fscores <- scores[colnames(sdata), ]

type <- c(rep("GENE", length(grep("ENSG", rownames(fulldata)))),
          rep("16S", length(grep("OTU", rownames(fulldata)))))

fdata <- data.frame(var = rownames(fulldata), type, dynamicColors, fscores)

(ts <- table(type, dynamicColors))
# Modules with bothg Genes or microorganisms
keep <- ts[1, ] != 0 & ts[2, ] != 0 & names(ts) != "grey"
out <- names(keep)[keep]
stopifnot(!is.null(out))

typecol <- rep("white", length(fdata$dynamicColors))
typecol[grep("OTU", rownames(fulldata))]<-"black"
typescorepc1<-rep("red", length(fdata$dynamicColors))
typescorepc1[which(fdata$V1<0)]<-"blue"
typescorepc2<-rep("red", length(fdata$dynamicColors))
typescorepc2[which(fdata$V2<0)]<-"blue"

triming <- function(num){
  for(i in 1:nrow(num)) {
    x <- num[i, ]
    trim = 0.05
    lo = quantile(x, trim)
    hi = quantile(x, 1 - trim)
    x[x < lo] = lo
    x[x > hi] = hi
    num[i, ] <- x
  }
  num
}

hmcol<-colorRampPalette(c("blue","white","red"))(552)

pdf("Figures/20160614_tTree_modules.pdf")
plotDendroAndColors(geneTree,
                    cbind(typecol,
                          dynamicColors,
                          typescorepc1,
                          typescorepc2), 
                    c("Omic","Module","F1","F2"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

labrow <- rep("", nrow(fulldata))

pdf("Figures/20190614_heatmap_2omics_with_modules_and_weights.pdf")
tf <- triming(fulldata)
heatmap.2(tf, col = hmcol, Rowv = as.dendrogram(geneTree), 
          density.info = "none", trace = "none", RowSideColors = dynamicColors,
          cexRow = 0.3, cexCol = 0.1, margins = c(5,10), labRow = labrow, 
          ColSideColors = labels2colors(meta_r[, 9]), key = FALSE)

# weights are proportional to increased in illeum !!

heatmap.2(tf, col = hmcol, Rowv = as.dendrogram(geneTree), 
          density.info = "none", trace = "none", RowSideColors = typecol,
          cexRow = 0.3, cexCol = 0.1, margins = c(5,10), labRow = labrow,
          ColSideColors = labels2colors(meta_r[, 9]), key = FALSE)

dev.off()

WriteXLS("fdata","Fdata.xls")
