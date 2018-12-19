library("ggforce")
library("RGCCA")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

# Select the Demographics
Demographics <- model_RGCCA(meta_r, c("ID","SEX", "Surgery", "Treatment", "Tobacco"))


# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "Demographics" = Demographics)
A <- clean_unvariable(A)

# The design
model <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0d <- subSymm(model, "16S", "Demographics", 1)
model0d <- subSymm(model0d, "RNAseq", "Demographics", 1)
model0di <- subSymm(model0d, 1, 1, 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

sgcca.centroid <- sgcca(
  A, C = model0d, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
sgcca.centroid$AVE$AVE_inner

saveRDS(sgcca.centroid, file = "model0_demographics.RDS")

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1],
  "Demographics" = sgcca.centroid$Y[["Demographics"]][, 1]
)


## Grouping of the variables ####
RNAseq1 <- samples$RNAseq
RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
microbiota2 <- sgcca.centroid$Y[["16S"]][, 2]
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


pdf(paste0("Figures/", today, "_RGCCA_plots_model0_demographics.pdf"))

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
result.out <- lapply(l, loo_model, model = model0d) 
saveRDS(result.out, "loo-model0d.RDS")
result.out <- lapply(l, loo_model, model = model0di) 
saveRDS(result.out, "loo-model0di.RDS")


# Bootstrap of sgcca
boot <- boot_sgcca(A, model0di, shrinkage, 1000)
saveRDS(boot, file = "bootstrap_model0di.RDS")
boot <- boot_sgcca(A, model0d, shrinkage, 1000)
saveRDS(boot, file = "bootstrap_model0d.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()
