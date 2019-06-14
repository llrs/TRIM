library("ggforce")
library("RGCCA")
library("BiocParallel")
library("integration")
library("fgsea")
library("broom")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("../intestinal_16S_RNAseq_integration/otus_table_norm_RNAseq.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r

Localization <- model_RGCCA(meta_r, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(meta_r, c("AgeDiag", "AGE_SAMPLE", "Transplant"))
Demographics <- model_RGCCA(meta_r, c("ID","SEX", "Surgery", "Treatment"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "Demographics" = Demographics,
          "Location" = Localization, "Time" = Time)

stopifnot(length(unique(vapply(A, nrow, numeric(1L)))) == 1)
A <- clean_unvariable(A)

saveRDS(A, "model3_TRIM.RDS")

# The design of model 3
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model3 <- subSymm(C, "16S", "Demographics", 1)
model3 <- subSymm(model3, "16S", "Time", 1)
model3 <- subSymm(model3, "16S", "Location", 1)
model3 <- subSymm(model3, "RNAseq", "Demographics", 1)
model3 <- subSymm(model3, "RNAseq", "Time", 1)
model3 <- subSymm(model3, "RNAseq", "Location", 1)
model3 <- subSymm(model3, "RNAseq", "16S", 1)
# C <- subSymm(C, "16S", "16S", 0)
# C <- subSymm(C, "RNAseq", "RNAseq", 0)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0, 1, 1, 1) # We guess a 0.1 for the RNAseq expression
# shrinkage[[2]] <- vapply(A[[2]], tau.estimate, numeric(1L))
# (min_shrinkage <- sapply(A, function(x) {
#   1 / sqrt(ncol(x))
# }))
# # Don't let the shrinkage go below the thershold allowed
# shrinkage <- ifelse(shrinkage < min_shrinkage | is.na(shrinkage), min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))
shrinkage[[2]] <- tau.estimate(A[[2]]) # 0.286506412433534

ncomp <- rep(2, length(A))


# A[1:2] <- lapply(A[1:2], function(x){scale2(x, bias = TRUE)})
sgcca.centroid <- sgcca(
  A, model3, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
saveRDS(sgcca.centroid, "sgcca_model3.RDS")
sgcca.centroid$AVE

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)

# designs <- weight_design(weights = 3, size = 5)
# keep <- check_design(designs)
# designs <- designs[keep]
# sgcca_custom <- function(x, ...) {
#   sgcca.centroid <- RGCCA::sgcca(
#     C = x,
#     scheme = "centroid",
#     scale = TRUE,
#     verbose = FALSE, ...)
#   sgcca.centroid$AVE[c("AVE_inner", "AVE_outer")]
# }
# ncomp <- rep(1, length(A))
# design_boot <- bplapply(designs, sgcca_custom, ncomp = ncomp,
#                         c1 = shrinkage, A = A, BPPARAM = bpparam())
# saveRDS(design_boot, "TRIM/intestinal_16S_RNAseq_metadb/sgcca_model3.RDS")

# Modify for a better usage
# w <- t(vapply(designs, function(x){x[upper.tri(x)]}, numeric(3L)))
# ind <- apply(which(upper.tri(designs[[1]]), arr.ind = TRUE), 1, 
#              paste0, collapse = "")
# colnames(w) <- paste0("var", ind)
# db <- t(vapply(design_boot, unlist, numeric(2L)))
# db2 <- cbind(db, w)
# db3 <- as.data.frame(db2)
# 
# lmM <- lm(AVE_inner~0+var12+var13+var23, data = db3)
# glance(lmM)
# tidy(lmM)

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1]
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


pdf(paste0("Figures/", today, "_RGCCA_plots_model3.pdf"))

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
comp1 <- sapply(sgcca.centroid$a, function(x) {  x[, 1]})
variables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {  x[, 2]})
variables_weight(comp2)

l <- looIndex(size(A))
loo_model <- loo_functions(A, shrinkage)
result.out <- lapply(l, loo_model, model = model3)
saveRDS(result.out, "loo-model3.RDS")

# Bootstrap of sgcca
boot <- boot_sgcca(A, model3, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model3.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()

