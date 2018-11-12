
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

# Tidify all the data of models:
names(model3_best$Y) <- names(model3$Y)
names(model3_besti$Y) <- names(model3$Y)
library("tidyr")
library("dplyr")

tidyer <- function(data, model, type) {
  d <- as.data.frame(data) %>% 
    as_tibble() %>% 
    mutate(Model = model) %>% 
    gather(Component, !!type, comp1:comp2)
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
m1GE <- merger(tidyer(model1$Y[[1]], "1", "GE"))
m1M <- merger(tidyer(model1$Y[[2]], "1", "M"))
m2GE <- merger(tidyer(model2$Y[[1]], "2", "GE"))
m2M <- merger(tidyer(model2$Y[[2]], "2", "M"))
m2bGE <- merger(tidyer(model2_best$Y[[1]], "2 best", "GE"))
m2bM <- merger(tidyer(model2_best$Y[[2]], "2 best", "M"))
m3GE <- merger(tidyer(model3$Y[[1]], "3", "GE"))
m3M <- merger(tidyer(model3$Y[[2]], "3", "M"))
m3bGE <- merger(tidyer(model3_best$Y[[1]], "3 best", "GE"))
m3bM <- merger(tidyer(model3_best$Y[[2]], "3 best", "M"))

inter <- intersect(colnames(m0GE), colnames(m0M))
inter <- grep("Rownames", inter, invert = TRUE, value = TRUE)

df <- rbind(
  merge(m0M, m0GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m1M, m1GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m2M, m2GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m2bM, m2bGE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m3M, m3GE, all.x = TRUE, all.y = TRUE, by = inter),
  merge(m3bM, m3bGE, all.x = TRUE, all.y = TRUE, by = inter)
)

library("ggplot2")
# Set theme without background on the labels
th <- theme_bw() +
  theme_update(strip.background = element_blank())
theme_set(th)
df <- as_tibble(df)
df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = as.factor(Sample))) +
  facet_wrap(~Model, scales = "free") + 
  guides(col = FALSE) +
  labs(title = "Samples by model", 
       subtitle = "Colored by sample",
       caption = "HSCT dataset")

# Check that the samples order doesn't change or something!! It doesn't look right
df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = IBD)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset")
df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = SESCD_local)) +
  facet_wrap(~Model, scales = "free") +
  labs(title = "Samples by model",
       caption = "HSCT dataset",
       col = "SESCD (local)")

df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = Exact_location)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Exact location")
df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = ID)) +
  facet_wrap(~Model, scales = "free") + 
  labs(title = "Samples by model",
       caption = "HSCT dataset", 
       col = "Patient ID")
df %>% 
  filter(Component == "comp1") %>% 
  ggplot() +
  geom_point(aes(GE, M, col = SEX)) +
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
       