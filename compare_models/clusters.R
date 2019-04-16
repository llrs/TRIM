library("integration")
library("gplots")
library("colorspace")
library("dfrtopics")
library("ggplot2")
library("ggridges")
library("org.Hs.eg.db")

# Evaluate all models ###
# model 0
folder0 <- "intestinal_16S_RNAseq_integration"
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))
model0_loo <- readRDS(file.path(folder0, "loo-model0.RDS"))
model0i <- readRDS(file.path(folder0, "sgcca_i.RDS"))
model0i_loo <- readRDS(file.path(folder0, "loo-model0i.RDS"))

# model 1 without interaction
folder1 <- "intestinal_16S_RNAseq_metadb"
model1 <- readRDS(file.path(folder1, "sgcca.RDS"))
model1_loo <- readRDS(file.path(folder1, "loo-model1.RDS"))
model1i <- readRDS(file.path(folder1, "sgccai.RDS"))
model1i_loo <- readRDS(file.path(folder1, "loo-model1i.RDS"))

# model 2 without interaction
model2 <- readRDS(file.path(folder1, "sgcca_model2.RDS"))
model2_loo <- readRDS(file.path(folder1, "loo-model2.RDS"))

# TODO create these files
# model2i <- readRDS(file.path(folder1, "sgcca_model2i.RDS"))
# model2i_loo <- readRDS(file.path(folder1, "loo-model2i.RDS"))

model2_best <- readRDS(file.path(folder1, "model2_best.RDS"))
model2_best_loo <- readRDS(file.path(folder1, "loo-model2_best.RDS"))

model2_besti <- readRDS(file.path(folder1, "model2_best_interaction.RDS"))
model2_besti_loo <- readRDS(file.path(folder1, "loo-model2_best_interaction.RDS"))

# model 3
model3 <- readRDS(file.path(folder1, "sgcca_model3.RDS"))
model3_loo <- readRDS(file.path(folder1, "loo-model3.RDS"))

model3_best <- readRDS(file.path(folder1, "model3_best.RDS"))
model3_best_loo <- readRDS(file.path(folder1, "loo-model3_best.RDS"))
model3_besti <- readRDS(file.path(folder1, "model3_best_interaction.RDS"))
model3_besti_loo <- readRDS(file.path(folder1, "loo-model3_best_interaction.RDS"))

model3_best2 <- readRDS(file.path(folder1, "model3_forced_interaction.RDS"))
model3_bestB <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))
model3_bestB_loo <- readRDS(file.path(folder1, "loo-model3_wo_forced_interaction.RDS"))


m0 <- weights_correlation(model0_loo)
m1.2 <- weights_correlation(model2_best_loo)
m2.1 <- weights_correlation(model3_best_loo)
m2.2 <- weights_correlation(model3_bestB_loo)


# Comes from models_sets
o <- read.csv("relevant_genes_different_models.csv")
o <- o[, 1]

mm0 <- m0[rownames(m0) %in% o, colnames(m0) %in% o]
mm2.1 <- m2.1[rownames(m2.1) %in% o, colnames(m2.1) %in% o]
mm2.2 <- m2.2[rownames(m2.2) %in% o, colnames(m2.2) %in% o]
mm1.2 <- m1.2[rownames(m1.2) %in% o, colnames(m1.2) %in% o]


mm0[is.na(mm0)] <- 0
mm2.1[is.na(mm2.1)] <- 0
mm2.2[is.na(mm2.2)] <- 0
mm1.2[is.na(mm1.2)] <- 0
h0 <- heatmap.2(mm0, scale = "none", trace = "none", main = "model 0",
          col = diverge_hsv(50))
h1.2 <- heatmap.2(mm1.2, scale = "none", trace = "none", main = "model 1.2",
          col = diverge_hsv(50))
h2.1 <- heatmap.2(mm2.1, scale = "none", trace = "none", main = "model 2.1",
          col = diverge_hsv(50))
h2.2 <- heatmap.2(mm2.2, scale = "none", trace = "none", main = "model 2.2",
          col = diverge_hsv(50))

#' Set the matrix of the heatmap as it is seeen to make it easier to see which
#' go with which
reorder_m <- function(x) {
  x$carpet[rev(rownames(x$carpet)), rev(colnames(x$carpet))]
}

split_groups <- function(m) {
  cols <- range(which(m[, 1] == 0))
  rows <- range(which(m[1, ] == 0))
  genes <- rownames(m)
  otus <- colnames(m)
  cluster_otus <- list(cluster1 = otus[seq_len(min(rows) - 1)], 
                       cluster2 = otus[seq(from = max(rows) + 1, to = ncol(m))])
  cluster_genes <- list(cluster1 = genes[seq_len(min(cols) - 1)], 
                       cluster2 = genes[seq(from = max(cols) + 1, to = nrow(m))])
  list(OTUs = cluster_otus, genes = cluster_genes)
}

m <- reorder_m(h0)
cluster_otus <- list(cluster1 = 1:45, cluster2 = 57:67)
cluster_genes <- list(cluster1 = 1:198, cluster2 = 261:372)
cluster_00 <- list(OTUs = cluster_otus, genes = cluster_genes)

cluster_0 <- split_groups(h0$carpet)
cluster_1.2 <- split_groups(h1.2$carpet)
m <- reorder_m(h2.2)
cluster_2.2 <- list(OTUs = list(cluster1 = 1:38, cluster2 = 39:67),
                    genes = list(cluster1 = , cluster2 = ))

(r_sign <-  cor_sign(158))

cor_df <- function(x, r_sign = NULL) {
  if (is.null(r_sign)) {
    r_sign <- 0
  }
  values <- which(abs(x) > r_sign, arr.ind = TRUE)
  
  d <- apply(values, 1, function(y) {
    r <- x[y[1], y[2]]
    c(r = r, pvalue = pvalue(r, 158))
  })
  values <- as.data.frame(values)
  values$row <- rownames(values)
  values$col <- colnames(x)[values$col]
  
  d <- t(d)
  cbind(values, d)
}

df0 <- cor_df(mm0, r_sign)
df1.2 <- cor_df(mm1.2, r_sign)
df2.2 <- cor_df(mm2.2, r_sign)

map_genes <- function(x) {
  mutate(x, Ensembl = trimVer(col),
                Symbol = mapIds(org.Hs.eg.db, keys = Ensembl, 
                                keytype = "ENSEMBL", 
                                column = "SYMBOL", multiVals = "first"),
                entrez = mapIds(org.Hs.eg.db, keys = Ensembl,
                                keytype = "ENSEMBL", 
                                column = "ENTREZID", multiVals = "first"))
}

df0 <- map_genes(df0)
df1.2 <- map_genes(df1.2)
df2.2 <- map_genes(df2.2)

write.csv(df1.2, file = "significant_LOO_correlations_model1.2.csv", na = "", 
          row.names = FALSE)
write.csv(df2.2, file = "significant_LOO_correlations_model2.2.csv", na = "", 
          row.names = FALSE)


m0_df <- gather_matrix(m0, col_names = c("OTUs", "Genes", "correlation"))
m0_df <- m0_df[!is.na(m0_df$correlation), ]
m1.2_df <- gather_matrix(m1.2, col_names = c("OTUs", "Genes", "correlation"))
m1.2_df <- m1.2_df[!is.na(m1.2_df$correlation), ]
m2.2_df <- gather_matrix(m2.2, col_names = c("OTUs", "Genes", "correlation"))
m2.2_df <- m2.2_df[!is.na(m2.2_df$correlation), ]

m0_df <- cbind(m0_df, model = "0")
m1.2_df <- cbind(m1.2_df, model = "1.2")
m2.2_df <- cbind(m2.2_df, model = "2.2")
m_df <- rbind(m0_df, m1.2_df, m2.2_df)
m_df <- mutate(m_df, keep = if_else(OTUs %in% o & Genes  %in% o, 1, 0))
ggplot(m_df) +
  geom_density(aes(x = correlation, y = ..scaled.., group = model, color = model)) +
  ggtitle("All")
ggplot(m_df[m_df$keep == 1, ]) +
  geom_density(aes(x = correlation, y = ..scaled.., group = model, color = model)) +
  ggtitle("subset")
ggplot(m_df) +
  geom_density_ridges(aes(correlation, y = model, color = model))
