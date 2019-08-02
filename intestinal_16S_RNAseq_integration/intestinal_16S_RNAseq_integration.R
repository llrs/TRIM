library("ggforce")
library("RGCCA")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")


# Save
otus_table_i <- readRDS("otus_table.RDS")
otus_n <- readRDS("otus_table_norm_RNAseq.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS( "meta.RDS")

# Prepare input for the sgcca function
A <- list("RNAseq" = t(expr), "16S" = t(otus_n))
A <- clean_unvariable(A)
saveRDS(A, file = "TRIM.RDS")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0 <- subSymm(C, "16S", "RNAseq", 1)
model0i <- subSymm(model0, 1, 1, 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0) # We guess a 0.1
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))
sgcca.centroid <- sgcca(
  A, C = model0i, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
beepr::beep()
saveRDS(sgcca.centroid, "sgccai.RDS")
sgcca.centroid <- sgcca(
  A, C = model0, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
saveRDS(sgcca.centroid, file = "sgcca.RDS")

# Find the direction of the correlation
samples <- data.frame(
  RNAseq = sgcca.centroid$Y[[1]][, 1],
  Micro = sgcca.centroid$Y[[2]][, 1]
)
lmt <- lm(Micro ~ RNAseq, data = samples)

if (lmt$coefficients[2] > 0) {
  d <- c(1, 1)
} else if (lmt$coefficients[2] < 0) {
  d <- c(-1, -1)
}


dist <- apply(samples, 1, dist2d, d = d)
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, meta_r, "dist" = dist)

# Plots ####
pdf(paste0("Figures/", today, "_RGCCA_plots.pdf"))

# Plot if the coherence between samples has a specific pattern
ggplot(samples) +
  geom_point(aes(Patient_ID, log10(dist), col = Involved_Healthy)) +
  facet_grid(~ Time)

hist(samples$dist)

plot_samples(samples, colors)

variables <- variables(sgcca)

plot_variables(variables)

# Plot for the same component the variables of each block
comp1 <- lapply(sgcca.centroid$a, function(x) {x[, 1]})
variables_weight(comp1)

# Second component
comp2 <- lapply(sgcca.centroid$a, function(x) {x[, 2]})
variables_weight(comp2)

# Validate #### 
l <- looIndex(size(A))
loo_model <- loo_functions(A, shrinkage)
# SGCCA of the selected model leaving one sample each time out of order.
result.out <- lapply(l, loo_model, model = model0) 
saveRDS(result.out, "loo-model0.RDS")
result.out <- lapply(l, loo_model, model = model0i) 
saveRDS(result.out, "loo-model0i.RDS")

# Bootstrap of sgcca
boot <- boot_sgcca(A, model0, shrinkage, 1000)
saveRDS(boot, file = "bootstrap.RDS")
booti <- boot_sgcca(A, model0i, shrinkage, 1000)
saveRDS(booti, file = "bootstrapi.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)
boot_evaluate(booti$STAB)

dev.off()
