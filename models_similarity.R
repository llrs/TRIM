library("BioCor")
library("ggplot2")
library("dplyr")
library("ggrepel")

# Evaluate all models ###
# model 0
folder0 <- "intestinal_16S_RNAseq_integration"
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))
# model 1 without interaction
folder1 <- "intestinal_16S_RNAseq_metadb"
model1 <- readRDS(file.path(folder1, "sgcca.RDS"))
# model 2 without interaction
model2 <- readRDS(file.path(folder1, "sgcca_model2.RDS"))
model2_best <- readRDS(file.path(folder1, "model2_best.RDS"))
# model 3
model3 <- readRDS(file.path(folder1, "sgcca_model3.RDS"))
model3_best <- readRDS(file.path(folder1, "model3_best.RDS"))


model3_best2 <- readRDS(file.path(folder1, "model3_forced_interaction.RDS"))
model3_bestB <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))

models <- list(model0, model1, model2, model2_best, model3, model3_best, 
               model3_best2, model3_bestB)
names(models) <- c("model0", "model1", "model2", "model2_best", "model3", 
                   "model3_best", "model3_best2", "model3_bestB")

genes <- sapply(models, function(x){
  rownames(x$a[[1]])[x$a[[1]][, 1] != 0]
}, simplify = FALSE)

micro <- sapply(models, function(x){
  rownames(x$a[[2]])[x$a[[2]][, 1] != 0]
})


genesSim <- mpathSim(names(models), BioCor:::inverseList(genes), method = NULL)
microSim <- mpathSim(names(models), BioCor:::inverseList(micro), method = NULL)

hc <- hclust(as.dist(1 - genesSim))
plot(hc, main = "Similarities between models", sub = "by genes")
hc2 <- hclust(as.dist(1 - microSim))
plot(hc2, main = "Similarities between models", sub = "by OTUs")

all_list <- c(BioCor:::inverseList(genes), BioCor:::inverseList(micro))
allSim <- mpathSim(names(models), all_list, method = NULL)

hc3 <- hclust(as.dist(1 - allSim))
plot(hc3, main = "Similarities between models", sub = "by OTUs and genes")

hc3b <- hclust(as.dist(1 - (genesSim + microSim)))
plot(hc3b, main = "Similarities between models", sub = "by OTUs and genes")



MDS_micro <- cmdscale(1 - microSim)
MDS_genes <- cmdscale(1 - genesSim)
MDS_micro <- as.data.frame(MDS_micro) %>% 
  tibble::rownames_to_column()
MDS_genes <- as.data.frame(MDS_genes) %>% 
  tibble::rownames_to_column()
colnames(MDS_genes) <- c("rowname", "PC1", "PC2")
colnames(MDS_micro) <- c("rowname", "PC1", "PC2")
MDS_genes <- cbind(MDS_genes, Type = "Genes")
MDS_micro <- cbind(MDS_micro, Type = "OTUs")
MDS <- rbind(MDS_genes, MDS_micro)
MDS %>% filter(Type == "Genes") %>% 
  ggplot() +
  geom_label_repel(aes(PC1, PC2, label = rowname)) +
  ggtitle("Genes PCA")

MDS %>% mutate(Dim = paste0("PC1_", Type))

MDS2 <- data.frame(rowname = MDS$rowname[1:8], 
           PC1_genes = MDS$PC1[1:8],
           PC1_micro = MDS$PC1[9:16])
ggplot(MDS2) +
  geom_label_repel(aes(PC1_genes, PC1_micro, label = rowname)) +
  ggtitle("Mixing PCAs")
