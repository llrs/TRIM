library("tidyr")
library("dplyr")
library("ggplot2")

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

getAVEs <- function(x, pos = 1L) { 
  x$AVE$AVE_inner[pos]
}

sem <- function(x){
  sd(x)/length(x)
}

plotAVEs <- function(model, loo) {
  aves <- vapply(loo, getAVEs, numeric(1L))
  hist(aves, xlim = c(0, 1), main = "model")
  abline(v = model$AVE$AVE_inner[1])
}

m_sem <- function(model, loo) {
  aves <- vapply(loo, getAVEs, numeric(1L))
  paste0(signif(model$AVE$AVE_inner[1], 3), 
        " (", signif(mean(aves), 3),
        " \u00b1 ", scales::scientific(sem(aves), 3), ")")
}

m_sem(model0, model0_loo)
m_sem(model0i, model0i_loo)
m_sem(model1, model1_loo)
m_sem(model1i, model1i_loo)
m_sem(model2, model2_loo)
m_sem(model2_best, model2_best_loo)
m_sem(model2_besti, model2_besti_loo)
m_sem(model3, model3_loo)
m_sem(model3_best, model3_best_loo)
m_sem(model3_besti, model3_besti_loo)

model0i$AVE$AVE_inner-model0$AVE$AVE_inner
model1i$AVE$AVE_inner-model1$AVE$AVE_inner
model2_besti$AVE$AVE_inner-model2_best$AVE$AVE_inner
model3_besti$AVE$AVE_inner-model3_best$AVE$AVE_inner

# Tidify all the data of models
names(model2_best$Y) <- names(model2$Y)
names(model3_best$Y) <- names(model3$Y)
names(model3_besti$Y) <- names(model3$Y)
names(model2_best$a) <- names(model2$Y)
names(model3_best$a) <- names(model3$Y)
names(model3_besti$a) <- names(model3$Y)

tidyer <- function(data, model, type) {
  if ("comp1" %in% colnames(data)){
    d <- as.data.frame(data) %>% 
      as_tibble() %>% 
      mutate(Model = model) %>% 
      gather(Component, !!type, comp1:comp2)
  } else {
    d <- as.data.frame(data) %>% 
      as_tibble() %>% 
      mutate(Model = model) %>% 
      gather(Component, !!type, 1:2)
  }
  d$Rownames <- rep(rownames(data), 2)
  d
    # mutate(Sample = seq_len(n()))
  # Samples name could be important!!
}  

meta <- readRDS("intestinal_16S_RNAseq_metadb/meta.RDS")

merger <- function(data) {
  if ("GE" %in% colnames(data)) {
    merge(data, meta, by.x = "Rownames", by.y = "Sample Name_RNA")
  } else {
    merge(data, meta, by.x = "Rownames", by.y = "Seq_code_uDNA")
  }
}

m0GE <- merger(tidyer(model0$Y[[1]], "0", "GE"))
m0M <- merger(tidyer(model0$Y[[2]], "0", "M"))
m0iGE <- merger(tidyer(model0i$Y[[1]], "0 i", "GE"))
m0iM <- merger(tidyer(model0i$Y[[2]], "0 i", "M"))
m1GE <- merger(tidyer(model1$Y[[1]], "1", "GE"))
m1M <- merger(tidyer(model1$Y[[2]], "1", "M"))
m1iGE <- merger(tidyer(model1i$Y[[1]], "1 i", "GE"))
m1iM <- merger(tidyer(model1i$Y[[2]], "1 i", "M"))
m2GE <- merger(tidyer(model2$Y[[1]], "2", "GE"))
m2M <- merger(tidyer(model2$Y[[2]], "2", "M"))
m2bGE <- merger(tidyer(model2_best$Y[[1]], "2 best", "GE"))
m2bM <- merger(tidyer(model2_best$Y[[2]], "2 best", "M"))
m2biGE <- merger(tidyer(model2_besti$Y[[1]], "2 best i", "GE"))
m2biM <- merger(tidyer(model2_besti$Y[[2]], "2 best i", "M"))
m3GE <- merger(tidyer(model3$Y[[1]], "3", "GE"))
m3M <- merger(tidyer(model3$Y[[2]], "3", "M"))
m3bGE <- merger(tidyer(model3_best$Y[[1]], "3 best", "GE"))
m3bM <- merger(tidyer(model3_best$Y[[2]], "3 best", "M"))
m3biGE <- merger(tidyer(model3_besti$Y[[1]], "3 best i", "GE"))
m3biM <- merger(tidyer(model3_besti$Y[[2]], "3 best i", "M"))

inter <- intersect(colnames(m0GE), colnames(m0M))
inter <- grep("Rownames", inter, invert = TRUE, value = TRUE)

df <- rbind(
  merge(m0M, m0GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m0iM, m0iGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m1M, m1GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m1iM, m1iGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m2M, m2GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m2bM, m2bGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m2biM, m2biGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m3M, m3GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m3bM, m3bGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m3biM, m3biGE, all.x = TRUE, all.y = TRUE, by = inter)
)

# Set theme without background on the labels
theme_set(theme_bw())
theme_update(strip.background = element_blank())
df <- as_tibble(df)
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = as.factor(Sample_Code_uDNA))) +
  facet_wrap(~Model, scales = "free") + 
  guides(col = FALSE) +
  labs(title = "Samples by model", 
       subtitle = "Colored by sample",
       caption = "HSCT dataset")

# Check that the samples order doesn't change or something!! It doesn't look right
df %>% 
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = IBD)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = SESCD_local)) +
  facet_wrap(~Model, scales = "free") +
  labs(title = "Samples by model",
       caption = "HSCT dataset",
       col = "SESCD (local)")

df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  mutate(Ileum = case_when(Exact_location == "ILEUM" ~ "Ileum", 
                           !is.na(Exact_location) ~ "Colon")) %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Ileum)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Location")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  mutate(Ileum = case_when(Exact_location == "ILEUM" ~ "Ileum", 
                           !is.na(Exact_location) ~ "Colon")) %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Ileum)) +
  stat_ellipse(aes(GE, M, col = Ileum, group = Ileum), type = "norm", show.legend = FALSE) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Location")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = ID)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Patient ID")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = IBD)) +
  stat_ellipse(aes(GE, M, col = IBD, group = IBD), type = "norm", show.legend = FALSE) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "IBD")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>%
  ggplot() +
  geom_point(aes(GE, M, col = Time)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Time")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Time)) +
  stat_ellipse(aes(GE, M, col = Time, group = Time), type = "norm", show.legend = FALSE) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Time")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Transplant)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Transplant")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Transplant)) +
  stat_ellipse(aes(GE, M, col = Transplant, group = Transplant), type = "norm", show.legend = FALSE) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Transplant")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = SEX)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Sex")
df %>% 
  filter(!grepl(" i", Model)) %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = SEX)) +
  stat_ellipse(aes(GE, M, col = SEX, group = SEX), type = "norm", show.legend = FALSE) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Sex")


# Plot one dimension vs the other
inter2 <- grep("Model|Sample_Code_uDNA", inter, value = TRUE, invert = TRUE)
sumDplyr <- function(x){
  y <- x[!is.na(x)]
  unique(y)
}
GEs <- df %>% 
  filter(!grepl(" i", Model)) %>% 
  spread(Component, GE) %>% 
  group_by(Model, Sample_Code_uDNA) %>% 
  summarise(GE1 = sum(comp1, na.rm = TRUE),
            GE2 =  sum(comp2, na.rm = TRUE),
            ID = unique(ID),
            SESCD_local = unique(SESCD_local),
            SEX = unique(SEX),
            IBD = unique(IBD)) %>% 
  ungroup()
Ms <- df %>% 
  filter(!grepl(" i", Model)) %>% 
  spread(Component, M) %>% 
  group_by(Model, Sample_Code_uDNA) %>% 
  summarise(M1 = sum(comp1, na.rm = TRUE),
            M2 = sum(comp2, na.rm = TRUE),
            ID = unique(ID),
            SESCD_local = unique(SESCD_local),
            SEX = unique(SEX),
            IBD = unique(IBD)
            ) %>% 
  ungroup()
GEs %>% 
  # filter(Model == "3 best") %>% 
  ggplot() +
  geom_point(aes(GE1, GE2, col = IBD)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model", 
       subtitle = "Gene expression dimensions",
       caption = "HSCT dataset")
GEs %>% 
  # filter(Model == "3 best") %>% 
  ggplot() +
  geom_point(aes(GE1, GE2, col = SEX)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model", 
       subtitle = "Gene expression dimensions",
       caption = "HSCT dataset")
       
Ms %>% 
  # filter(Model == "3 best") %>% 
  ggplot() +
  geom_point(aes(M1, M2, col = ID)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       subtitle = "Microbiome dimensions",
       caption = "HSCT dataset")
Ms %>% 
  # filter(Model == "3 best") %>% 
  ggplot() +
  geom_point(aes(M1, M2, col = IBD)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       subtitle = "Microbiome dimensions",
       caption = "HSCT dataset")
       
Ms %>% 
  # filter(Model == "3 best") %>% 
  ggplot() +
  geom_point(aes(M1, M2, col = SEX)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       subtitle = "Microbiome dimensions",
       caption = "HSCT dataset")

# Look the weights ####
a0GE <- tidyer(model0$a[[1]], "0", "GE")
a0M <- tidyer(model0$a[[2]], "0", "M")
a0iGE <- tidyer(model0i$a[[1]], "0 i", "GE")
a0iM <- tidyer(model0i$a[[2]], "0 i", "M")
a1GE <- tidyer(model1$a[[1]], "1", "GE")
a1M <- tidyer(model1$a[[2]], "1", "M")
a1iGE <- tidyer(model1i$a[[1]], "1 i", "GE")
a1iM <- tidyer(model1i$a[[2]], "1 i", "M")
a2GE <- tidyer(model2$a[[1]], "2", "GE")
a2M <- tidyer(model2$a[[2]], "2", "M")
a2bGE <- tidyer(model2_best$a[[1]], "2 best", "GE")
a2bM <- tidyer(model2_best$a[[2]], "2 best", "M")
a2biGE <- tidyer(model2_besti$a[[1]], "2 best i", "GE")
a2biM <- tidyer(model2_besti$a[[2]], "2 best i", "M")
a3GE <- tidyer(model3$a[[1]], "3", "GE")
a3M <- tidyer(model3$a[[2]], "3", "M")
a3bGE <- tidyer(model3_best$a[[1]], "3 best", "GE")
a3bM <- tidyer(model3_best$a[[2]], "3 best", "M")
a3biGE <- tidyer(model3_besti$a[[1]], "3 best i", "GE")
a3biM <- tidyer(model3_besti$a[[2]], "3 best i", "M")

dfGE <- rbind(a0GE, a0iGE, a1GE, a1iGE, a2GE, a2bGE, a2biGE, a3GE, a3bGE, a3biGE)
dfM <- rbind(a0M, a0iM, a1M, a1iM, a2M, a2bM, a2biM, a3M, a3bM, a3biM)
keepGE <- dfGE %>% 
  # filter(!grepl(" i", Model)) %>% 
  filter(Component == "V1" & GE != 0) %>% 
  mutate(Presence = if_else(GE != 0, 1, 0)) %>% 
  select(-Component, -GE, Rownames) %>% 
  group_by(Rownames) %>% 
  spread(Model, Presence) %>% 
  as.data.frame()
rownames(keepGE) <- keepGE$Rownames
keepGE <- keepGE[, -grep("Rownames", colnames(keepGE))]
keepGE[is.na(keepGE)] <- 0

keepM <- dfM %>% 
  # filter(!grepl(" i", Model)) %>%
  filter(Component == "V1" & M != 0) %>% 
  mutate(Presence = if_else(M != 0, 1, 0)) %>% 
  select(-Component, -M, Rownames) %>% 
  group_by(Rownames) %>% 
  spread(Model, Presence) %>% 
  ungroup() %>% 
  as.data.frame()
rownames(keepM) <- keepM$Rownames
keepM <- keepM[, -grep("Rownames", colnames(keepM))]
keepM[is.na(keepM)] <- 0

library("UpSetR")
library("grid")
text_sizes <- c(1.3, 1.3, 1, 1, 1.5, 1.5)
dfGE %>% 
  filter(Component == "V1" & GE != 0) %>% 
  ggplot() +
  geom_density(aes(GE)) +
  facet_wrap("Model") +
  labs(title = "Distribution of the weights", xlab = "weights", 
       subtitle = "Gene expression")
dfM %>% 
  filter(Component == "V1" & M != 0) %>% 
  ggplot() +
  geom_density(aes(M)) +
  facet_wrap("Model") +
  labs(title = "Distribution of the weights", xlab = "weights", 
       subtitle = "Microbiome")

## Upset plots ####
upset(keepGE, order.by = "freq", nsets = 6,
      sets = rev(colnames(keepGE)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes, scale.sets = "identity")
grid.text("Genes shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
upset(keepM, order.by = "freq", nsets = 6, 
      sets = rev(colnames(keepM)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes)
grid.text("OTUs shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))

# Testing genes in model 1 to 3 best 
se <- apply(keepGE, 1, function(x){all(x[2:6] != 0 & x[1] == 0)})
library("org.Hs.eg.db")
g <- mapIds(org.Hs.eg.db, keys = integration::trimVer(names(se)[se]), column = "SYMBOL", keytype = "ENSEMBL")
g <- unique(g[!is.na(g)]) # Calco related S100A6
int_genes <- "^S100A|^NOD|DUOX|^CA[1:9]+|^CEACAM|^REG"
grep(int_genes, g, value = TRUE) # S100A6

# 0 to 3
se <- apply(keepGE, 1, function(x){all(x[1:5] != 0 & x[6] == 0)})
g2 <- mapIds(org.Hs.eg.db, keys = integration::trimVer(names(se)[se]), column = "SYMBOL", keytype = "ENSEMBL")
g2 <- unique(g2[!is.na(g2)]) #
grep(int_genes, g2, value = TRUE) # DUOXA[12] DUOX2 S100A11 S100A13

# 0 to 2 best
se <- apply(keepGE, 1, function(x){all(x[1:4] != 0 & x[5:6] == 0)})
g3 <- mapIds(org.Hs.eg.db, keys = integration::trimVer(names(se)[se]), column = "SYMBOL", keytype = "ENSEMBL")
g3 <- unique(g3[!is.na(g3)]) 
grep(int_genes, g3, value = TRUE) # S100A8 S100A9 S100A12

# Just 3 best
se <- apply(keepGE, 1, function(x){all(x[6] != 0 & x[1:5] == 0)})
g4 <- mapIds(org.Hs.eg.db, keys = integration::trimVer(names(se)[se]), column = "SYMBOL", keytype = "ENSEMBL")
g4 <- unique(g4[!is.na(g4)]) # DUOX1
grep(int_genes, g4, value = TRUE) # Kevin genes

g0 <- trimVer(rownames(keepGE))
g0 <- mapIds(org.Hs.eg.db, keys = g0, column = "SYMBOL", keytype = "ENSEMBL")
g0 <- unique(g0[!is.na(g0)])

# Looking for functional annotation 
entrezID <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(keepGE)), keytype = "ENSEMBL", 
                   column = "ENTREZID")

paths2genes <- integration:::access_reactome()
genes <- unlist(paths2genes, use.names = FALSE)
pathways <- rep(names(paths2genes), lengths(paths2genes))
desc <- mapIds(reactome.db, keys = pathways, column = "PATHNAME", keytype = "PATHID")
desc <- gsub("Homo sapiens: ", "", desc)
T2D <- unique(cbind.data.frame(pathways, desc))
T2G <- cbind.data.frame(pathways, genes)
library("clusterProfiler")
enrichment_models <- lapply(keepGE, function(x){
  y <- entrezID[x == 1]
  enrich <- enricher(gene = y, TERM2GENE = T2G, TERM2NAME = T2D)
  as.data.frame(enrich)
})

folder <- "Summary_out"
sapply(names(enrichment_models), function(x) {
  name <- paste0("enrichment_GE_model_", gsub(" ", "_", x), ".csv")
  write.csv(enrichment_models[[x]], 
            file = file.path(folder, name), 
            row.names = FALSE)
})
a <- sapply(enrichment_models, function(x) {
  x[, 1]
})
# Upset of the enrichment ####
paths <- unique(unlist(a, use.names = FALSE))
b <- sapply(a, function(x){
  as.numeric(paths  %in% x)
})
rownames(b) <- paths
upset(as.data.frame(b), order.by = "freq", nsets = 6, 
      sets = rev(colnames(b)), keep.order = TRUE,
      line.size = NA, text.scale = text_sizes)
grid.text("Pathways shared in models", x = 0.65, y = 0.95, gp = gpar(fontsize = 20))
m <- list(model0, model1, model2, model2_best, model3, model3_best)
names(m) <- colnames(keepGE)
sapply(names(m), function(x) {
  name <- paste0("weights_GE_model_", gsub(" ", "_", x), ".csv")
  write.csv(integration::weights(m[[x]]), 
            file = file.path(folder, name), 
            row.names = FALSE)
})
taxa <- readRDS("intestinal_16S_RNAseq_metadb/otus_tax.RDS")
sapply(names(m), function(x) {
  name <- paste0("weights_M_model_", gsub(" ", "_", x), ".csv")
  write.csv(integration::weights_otus(m[[x]], taxa), 
            file = file.path(folder, name), 
            row.names = FALSE, na = "")
})
# TODO look taxa of otus
tax_i <- readRDS("intestinal_16S_RNAseq_metadb/otus_tax.RDS")
# TODO use BioCor to compare subgroups
library("BioCor")
out <- dfGE %>% 
  filter(GE != 0, Component == "V1") %>% 
  mutate(Genes = mapIds(org.Hs.eg.db, keys = trimVer(Rownames), keytype = "ENSEMBL", column = "ENTREZID")) %>% 
  as.data.frame()
model2gene <- split(out$Genes, out$Model)
model2gene <- lapply(model2gene, unique)
gsc <- GSEAdv::as.GeneSetCollection(GSEAdv:::inverseList(paths2genes))
o <- mclusterGeneSim(clusters = model2gene, info = gsc, method = c("max", "BMA"))
saveRDS(model2gene, file = "models2genes.RDS")

specificGenes <- lapply(seq_len(6), function(x){
  specific <- keepGE[, x] != 0
  others <- apply(keepGE[, seq_len(6)[-x]], 1, function(z){
    all(z == 0)
  })
  y <- trimVer(rownames(keepGE)[specific & others])
  y2 <- mapIds(org.Hs.eg.db, keys = y, keytype = "ENSEMBL", column = "ENTREZID")
  y2 <- unique(y2)
  y2[!is.na(y2)]
})
names(specificGenes) <- colnames(keepGE)
o2 <- mclusterGeneSim(clusters = specificGenes[1:4], info = gsc, method = c("max", "BMA"))
o3 <- mclusterGeneSim(clusters = specificGenes, info = gsc, method = c("max", "BMA"))
o3 <- readRDS("compareClusters.RDS")
o3
o2 <- mclusterSim(clusters = specificGenes[1:4], info = gsc, method = "max")
o3 <- mclusterSim(clusters = specificGenes, info = gsc, method = "max")
save(specificGenes, gsc, file = "specificGenes_gsc_BioCor.RData")


specificOTUs <- vapply(seq_len(6), function(x){
  specific <- keepM[, x] != 0
  others <- apply(keepM[, seq_len(6)[-x]], 1, function(z){
    all(z == 0)
  })
  sum(specific & others)
}, numeric(1L))
names(specificOTUs) <- colnames(keepM)
