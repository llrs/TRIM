library("ggforce")
library("RGCCA")
library("BiocParallel")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("../intestinal_16S_RNAseq_integration/otus_table_norm_RNAseq.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r
keepCol <- sapply(metadb, is.factor)
nam <- c(
  "Exact_location", # Segment of the sample
  # superseeded by SESCD 
  # "Active_area", # Health stage of the sample
  # "IBD", # Disease or control
  "AGE_SAMPLE", # Age
  # "diagTime", # Time with disease
  # Not really needed induced by diagTime and age sample
  "AgeDiag", # Age at which the disease was diagnositcated
  "Transplant", # Stage of the treatment
  "ID", # Patient
  # "SESCD_local", # Clinical score
  "Treatment", # Further complications
  "Surgery", # Up to surgery?
  "SEX" # Male/female
) 

metadb <- model_RGCCA(meta_r, nam) 

# Set metadb with a sigle variable with several options
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "metadata" = metadb)
A <- clean_unvariable(A)

# The design
model <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model1 <- subSymm(model, "16S", "metadata", 1)
model1 <- subSymm(model1, "RNAseq", "metadata", 1)
model2 <- subSymm(model1, "RNAseq", "16S", 1)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]]) # 0.286506412433534
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

sgcca.centroid <- sgcca(
  A, model2, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
sgcca.centroid$AVE

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)
saveRDS(sgcca.centroid, file = "sgcca_model2.RDS")

designs <- weight_design(weights = 11, size = length(A))
keep <- check_design(designs)

designs <- designs[keep]
sgcca_custom <- function(x, ...) {
  sgcca.centroid <- RGCCA::sgcca(
    C = x,
    scheme = "centroid",
    scale = FALSE,
    verbose = FALSE, ...)
  sgcca.centroid$AVE[c("AVE_inner", "AVE_outer")]
}

Ab <- lapply(A, scale2)
# design_boot <- bplapply(designs, sgcca_custom, ncomp = ncomp, 
# shrinkage = shrinkage, A = A, BPPARAM = bpparam())
design_boot <- lapply(designs, sgcca_custom, c1 = shrinkage, A = Ab)
# Modify for a better usage
w <- t(vapply(designs, function(x){x[upper.tri(x)]}, numeric(3L)))
ind <- apply(which(upper.tri(designs[[1]]), arr.ind = TRUE), 1, 
             paste0, collapse = "")
colnames(w) <- paste0("var", ind)
db <- t(vapply(design_boot, unlist, numeric(2L)))
db2 <- cbind(db, w)
db3 <- as.data.frame(db2)
saveRDS(db3, "designs_boot_model2.RDS")
design_boot <- readRDS("designs_boot_model2.RDS")

model2b <- design_boot[which.max(design_boot$AVE_inner), 3:5]
model2b <- symm(model1, model2b)
# Continue with the normal model 2
samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1],
  "metadata" = sgcca.centroid$Y[["metadata"]][, 1]
)

## Grouping of the variables ####
RNAseq1 <- samples$RNAseq
RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
microbiota2 <- sgcca.centroid$Y[["16S"]][, 2]
microbiota1 <- samples$Micro

names(RNAseq1) <- rownames(samples)
names(microbiota1) <- rownames(samples)
names(RNAseq2) <- rownames(samples)
names(microbiota2) <- rownames(samples)
groups <- split(rownames(samples), as.factor(meta_r$HSCT_responder))
# First dimension seems to capture well the
fgsea(groups, RNAseq1, nperm = 1000)
fgsea(groups, microbiota1, nperm = 1000)
# Further dimensions
fgsea(groups, RNAseq2, nperm = 1000)
fgsea(groups, microbiota2, nperm = 1000)


pdf(paste0("Figures/", today, "_RGCCA_plots_model2.pdf"))

km <- kmeans(samples[, c("RNAseq", "Micro")], 2, nstart = 2)
plot(samples[, c("RNAseq", "Micro")], col = km$cluster, 
     main = "K-clustering (2 groups)")

## Plotting ####
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, droplevels(meta_r))
plot_samples(samples, colors)

variables <- variables(sgcca.centroid)

plot_variables(variables)

# Plot PCAs
plot_interesting(variables, meta_r, expr, otus_table_i)

# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x) {x[, 1]})
variables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {x[, 2]})
variables_weight(comp2)

# Validate ####

model2_best <- matrix(0, nrow = 3, ncol = 3)
model2_best <- symm(model2_best, unlist(db3[which.max(db3$AVE_inner), 3:5]))
colnames(model2_best) <- colnames(model2)
rownames(model2_best) <- rownames(model2)

model2_best_sgcca <- sgcca(A, C = model2b, verbose = FALSE, c1 = shrinkage[1:3], ncomp = ncomp[1:3])
saveRDS(model2_best_sgcca, "model2_best.RDS")

model2_besti <- subSymm(model2b, 1, 1, 1)
model2_best_interaction_sgcca <- sgcca(A, C = model2_besti, verbose = FALSE, c1 = shrinkage[1:3], ncomp = ncomp[1:3])
saveRDS(model2_best_interaction_sgcca, "model2_best_interaction.RDS")

l <- looIndex(size(A))
loo_model <- loo_functions(A, shrinkage)
result.out <- lapply(l, loo_model, model = model2_best)
saveRDS(result.out, "loo-model2_best.RDS")
result.out <- lapply(l, loo_model, model = model2_besti)
saveRDS(result.out, "loo-model2_best_interaction.RDS")

result.out <- lapply(l, loo_model, model = model2)
saveRDS(result.out, "loo-model2.RDS")

# Bootstrap of sgcca
boot <- boot_sgcca(A, model2, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model1.2.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()

