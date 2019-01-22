library("ggforce")
library("RGCCA2")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
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
  "SEX"#, # Male/female
  # "Tobacco"
) 

metadb_0 <- model_RGCCA(meta_r, nam) 
metadb_0$AgeDiag[is.na(metadb_0$AgeDiag)] <- metadb_0$AGE_SAMPLE[is.na(metadb_0$AgeDiag)]

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "metadata" = metadb_0)
A <- clean_unvariable(A)
saveRDS(A, file = "TRIM.RDS")

# The design
model <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model1 <- subSymm(model, "16S", "metadata", 1)
model1 <- subSymm(model1, "RNAseq", "metadata", 1)
model1i <- subSymm(model1, 1, 1, 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0, 1) # We guess a 0.249488046688595 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(1, length(A))
sgcca.centroid <- sgcca(
  A, C = model1i, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
saveRDS(sgcca.centroid, file = "sgccai.RDS")
sgcca.centroid <- sgcca(
  A, C = model1, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
sgcca.centroid$AVE$AVE_inner
# 
# sgcca.factorial <- sgcca(
#   A, C, c1 = shrinkage,
#   ncomp = ncomp,
#   scheme = "factorial",
#   scale = TRUE,
#   verbose = FALSE
# )
# sgcca.factorial <- improve.sgcca(sgcca.factorial, names(A))
# 
# sgcca.horst <- sgcca(
#   A, C, c1 = shrinkage,
#   ncomp = ncomp,
#   scheme = "horst",
#   scale = TRUE,
#   verbose = FALSE
# )
# sgcca.horst <- improve.sgcca(sgcca.horst, names(A))

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)
saveRDS(sgcca.centroid, file = "sgcca.RDS")

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1],
  "metadata" = sgcca.centroid$Y[["metadata"]][, 1]
)


## Grouping of the variables ####
RNAseq1 <- samples$RNAseq
RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
microbiota2 <- sgcca.centroid$Y[["Micro"]][, 2]
microbiota1 <- samples$microbiota

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


pdf(paste0("Figures/", today, "_RGCCA_plots.pdf"))

km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
plot(samples[, c("RNAseq", "microbiota")], col = km$cluster, 
     main = "K-clustering (2 groups)")

## Plotting ####
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, droplevels(meta_r))
plot_samples(samples, colors)

variables <- variables(sgcca.centroid)

plot_variables(variables)

# Plot PCAs
plot_interesting(plot_interesting, meta_r, expr, otus_table_i)

# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x) {x[, 1]})
variables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {x[, 2]})
variables_weight(comp2)


# Validate ####

l <- looIndex(size(A))
loo_model <- loo_functions(A, shrinkage)
# SGCCA of the selected model leaving one sample each time out of order.
result.out <- lapply(l, loo_model, model = model1) 
saveRDS(result.out, "loo-model1.RDS")
result.out <- lapply(l, loo_model, model = model1i) 
saveRDS(result.out, "loo-model1i.RDS")


# Bootstrap of sgcca
boot <- boot_sgcca(A, model1, shrinkage, 1000)
saveRDS(boot, file = "bootstrap.RDS")
booti <- boot_sgcca(A, model1i, shrinkage, 1000)
saveRDS(booti, file = "bootstrapi.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)
boot_evaluate(booti$STAB)

dev.off()

