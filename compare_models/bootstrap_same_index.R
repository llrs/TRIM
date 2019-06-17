library("patchwork")
library("bayesboot")
library("ggplot2")
library("RGCCA")
library("integration")

# Models ###
folder0 <- "../intestinal_16S_RNAseq_integration"
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))

# model 1 without interaction
folder1 <- "../intestinal_16S_RNAseq_metadb"
# model 2 without interaction
model1.2 <- readRDS(file.path(folder1, "model2_best.RDS"))
# model 3
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))

# Read initial data ####
RNAseq <- readRDS(file.path(folder1, "expr.RDS"))
microbiome <- readRDS("../intestinal_16S_RNAseq_integration/otus_table_norm_RNAseq.RDS")
meta <-  readRDS(file.path(folder1, "meta.RDS"))
A <- list("RNAseq" = t(RNAseq), "16S" = t(microbiome), "metadata" = meta)

# Build the data ####
nam <- c("Exact_location", "AGE_SAMPLE", "AgeDiag", "Transplant", "ID",
         "Treatment", "Surgery", "SEX")
A[[3]]$AgeDiag[is.na(A[[3]]$AgeDiag)] <- 0 # Time has NA values
A1.2 <- A
A1.2[[3]] <- model_RGCCA(A1.2[[3]], nam)
A1.2 <- clean_unvariable(A1.2)

A2.2 <- A

Localization <- model_RGCCA(A2.2[[3]], c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(A2.2[[3]], c("AgeDiag", "AGE_SAMPLE", "Transplant"))
Demographics <- model_RGCCA(A2.2[[3]], c("ID","SEX", "Surgery", "Treatment"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values

A2.2[[3]] <- Demographics
A2.2[[4]] <- Localization
A2.2[[5]] <- Time
A2.2 <- clean_unvariable(A2.2)
names(A2.2) <- c("RNAseq", "16S", "Demographics", "Location", "Time")

# Boots ####
seed <- 487178
set.seed(seed)
index <- vector("list", length = 1000)
for (i in seq_len(1000)) {
  index[[i]] <- sample(nrow(A[[1]]), replace = TRUE)
}
saveRDS(index, "index_locale.RDS")
index <- readRDS("index_locale.RDS")

base_boot <- function(index, A, C) {
  STAB <- vector("list", length = length(A))
  AVE <- vector("numeric", length = 2)
  names(AVE) <- c("inner", "outer")
  names(STAB) <- names(A)
  
  A <- subsetData(A, index)
  
  try({
    res <- sgcca(A, C, scheme = "centroid", scale = TRUE)
    AVE["inner"] <- res$AVE$AVE_inner
    AVE["outer"] <- res$AVE$AVE_outer
    for (j in seq_along(A)) {
      STAB[[j]] <- res$a[[j]][, 1]
    }
    
  }, silent = TRUE)
  list(AVE = AVE, STAB = STAB)
}

boot0 <- lapply(index, base_boot, A = A[1:2], C = model0$C)
saveRDS(boot0, "boot0.RDS")
boot1.2 <- lapply(index, base_boot, A = A1.2, C = model1.2$C)
saveRDS(boot1.2, "boot1.2.RDS")
boot2.2 <- lapply(index, base_boot, A = A2.2, C = model2.2$C)
saveRDS(boot2.2, "boot2.2.RDS")