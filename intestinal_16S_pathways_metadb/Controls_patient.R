library("ggforce")
library("RGCCA")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

IBD <- meta_r$IBD == "CD"
meta_r <- meta_r[!IBD, ]
expr <- expr[, !IBD]
otus_table_i <- otus_table_i[, !IBD]

gsv <- readRDS("GSV.RDS")

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r
keepCol <- sapply(metadb, is.factor)
nam <- c(
  "Exact_location", # Segment of the sample
  # superseeded by SESCD 
  # "Active_area", # Health stage of the sample
  # "IBD", # Disease or control
  "AGE_SAMPLE", # Age
  "diagTime", # Time with disease
  # Not really needed induced by diagTime and age sample
  "AgeDiag", # Age at which the disease was diagnositcated
  "Transplant", # Stage of the treatment
  "ID", # Patient
  "SESCD_local", # Clinical score
  "Treatment", # Further complications
  "Surgery", # Up to surgery?
  "SEX" # Male/female
) 
metadb <- model_RGCCA(meta_r, nam) 
# Set metadb with a sigle variable with several options
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(gsv), "16S" = t(otus_table_i), "metadata" = metadb)
A <- clean_unvariable(A)

saveRDS(A, file = "TRIM_Controls.RDS")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "metadata", 1)
C <- subSymm(C, "RNAseq", "metadata", 1)


# We cannnot comput eht tau.estimate for A[[1]]
(shrinkage <- sapply(A, tau.estimate))
shrinkage["metadata"] <- 1
# shrinkage <- c(0.25670333, 0, 1) # We guess a 0.1 for the RNAseq expression
# shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

sgcca.centroid <- sgcca(
  A, C, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
saveRDS(sgcca.centroid, file = "Controls.RDS")

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1],
  "metadata" = sgcca.centroid$Y[["metadata"]][, 1]
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

km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
pdf(paste0("Figures/", today, "_RGCCA_plots_controls.pdf"))
plot(samples[, c("RNAseq", "microbiota")], col = km$cluster)

## Plotting results ####
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

# Bootstrap of sgcca
boot <- boot_sgcca(A, C, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_Controls.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()
