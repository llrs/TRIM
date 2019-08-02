library("RGCCA")
library("integration")
library("WGCNA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load
expr <- readRDS("expr.RDS")

pseudo_otus <- expr[sample(nrow(expr), 500), ]
# pseudo_otus <- apply(pseudo_otus, 1, sample, size = 158)
pseudo_otus <- apply(pseudo_otus, 1, function(x){
  runif(length(x)) + x
})
colnames(pseudo_otus) <- paste0("OTU_", seq_len(ncol(pseudo_otus)))
sd2 <- apply(pseudo_otus, 2, range)
sd0 <- apply(t(expr), 2, range)
A <- list("RNAseq" = t(expr), "16S" = pseudo_otus)

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0 <- subSymm(C, "16S", "RNAseq", 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0) # We guess a 0.1
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))


ncomp <- rep(2, length(A))
sgcca.centroid <- sgcca(
  A, C = model0, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)

sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
sgcca.centroid$AVE
plot(sgcca.centroid$Y[[1]][, 1], sgcca.centroid$Y[[2]][, 1])
summary(sgcca.centroid$a[[1]][, 1] != 0)

wgcna <- cbind(A[[1]], A[[2]])
selected <- sapply(sgcca.centroid$a, function(x){
  rownames(x)[x[, 1] != 0]
})
wgnca_small <- wgcna[, unlist(selected)]
# WGCNA ####

allowWGCNAThreads(3)
input <- wgnca_small
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
adjacency <- adjacency(input, power = 15)
dissTOM <- 1 - TOMsimilarity(adjacency)
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 10
dynamicMods <- cutreeDynamic(dendro = geneTree,  method = "tree", 
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

type <- c(rep("GENE", length(grep("ENSG", rownames(fulldata)))),
          rep("16S", length(grep("OTU", rownames(fulldata)))))

fdata <- data.frame(var = rownames(fulldata), type, dynamicColors)

(ts <- table(type, dynamicColors))
# Modules with bothg Genes or microorganisms
keep <- ts[1, ] != 0 & ts[2, ] != 0 & colnames(ts) != "grey"
names(keep)[keep]


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

pdf("Figures/Tree_modules.pdf")
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

