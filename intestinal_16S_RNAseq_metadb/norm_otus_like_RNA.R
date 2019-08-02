library("ggforce")
library("RGCCA")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("../intestinal_16S_RNAseq_integration/otus_table_norm_RNAseq.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
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
shrinkage[[2]] <- tau.estimate(A[[2]]) # from 0.286506412433534 to 0.306689360110999

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

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "Micro" = sgcca.centroid$Y[["16S"]][, 1]
)


## Grouping of the variables ####
RNAseq1 <- samples$RNAseq
RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
microbiota2 <- sgcca.centroid$Y[["Micro"]][, 2]
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



km <- kmeans(samples[, c("RNAseq", "Micro")], 2, nstart = 2)
plot(samples[, c("RNAseq", "Micro")], col = km$cluster, 
     main = "K-clustering (2 groups)")
plot(samples[, c("RNAseq", "Micro")], 
     col = as.factor(ifelse(meta_r$Exact_location == "ILEUM", "ILEUM", "COLON")),
     pch = as.numeric(as.factor(meta_r$IBD)) + 16)

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

