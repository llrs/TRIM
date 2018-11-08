
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
