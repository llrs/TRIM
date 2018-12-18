library("ggforce")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r

Localization <- model_RGCCA(meta_r, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(meta_r, c("AgeDiag", "AGE_SAMPLE", "Transplant"))
Demographics <- model_RGCCA(meta_r, c("ID","SEX", "Surgery", "Treatment"))

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "Demographics" = Demographics,
          "Location" = Localization, "Time" = Time)

stopifnot(length(unique(vapply(A, nrow, numeric(1L)))) == 1)
A <- lapply(A, function(x) {
  k <- unlist(apply(x, 2, sd, na.rm = FALSE))
  x[, k != 0]}
)

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
shrinkage <- c(0.25670333, 0, 1, 1, 1) # We guess a 0.1 for the RNAseq expression
# shrinkage[[2]] <- vapply(A[[2]], tau.estimate, numeric(1L))
# (min_shrinkage <- sapply(A, function(x) {
#   1 / sqrt(ncol(x))
# }))
# # Don't let the shrinkage go below the thershold allowed
# shrinkage <- ifelse(shrinkage < min_shrinkage | is.na(shrinkage), min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))
shrinkage[[2]] <- tau.estimate(A[[2]])

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
saveRDS(sgcca.centroid, file = "sgcca_model3.RDS")

designs <- weight_design(weights = 3, size = 5)
keep <- check_design(designs)
library("BiocParallel")
designs <- designs[keep]
sgcca_custom <- function(x, ...) {
  sgcca.centroid <- RGCCA::sgcca(
    C = x,
    scheme = "centroid",
    scale = TRUE,
    verbose = FALSE, ...)
  sgcca.centroid$AVE[c("AVE_inner", "AVE_outer")]
}
ncomp <- rep(1, length(A))
# design_boot <- bplapply(designs, sgcca_custom, ncomp = ncomp,
#                         c1 = shrinkage, A = A, BPPARAM = bpparam())
# saveRDS(design_boot, "TRIM/intestinal_16S_RNAseq_metadb/sgcca_model3.RDS")

# Modify for a better usage
w <- t(vapply(designs, function(x){x[upper.tri(x)]}, numeric(3L)))
ind <- apply(which(upper.tri(designs[[1]]), arr.ind = TRUE), 1, 
             paste0, collapse = "")
colnames(w) <- paste0("var", ind)
db <- t(vapply(design_boot, unlist, numeric(2L)))
db2 <- cbind(db, w)
db3 <- as.data.frame(db2)
library("broom")
lmM <- lm(AVE_inner~0+var12+var13+var23, data = db3)
glance(lmM)
tidy(lmM)

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "microbiota" = sgcca.centroid$Y[["16S"]][, 1]
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


pdf(paste0("Figures/", today, "_RGCCA_plots_model3.pdf"))

km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
plot(samples[, c("RNAseq", "microbiota")], col = km$cluster, 
     main = "K-clustering (2 groups)")

## Plotting ####
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, droplevels(meta_r))
samples$Patient_ID <- as.factor(samples$Patient_ID)
samples$Sample_Code <- as.character(samples$Sample_Code)

# Labels of the samples
label <- strsplit(as.character(samples$`Sample Name_RNA`), split = "-")
labels <- sapply(label, function(x) {
  if (length(x) == 5) {
    x[5]
  }
  else if (length(x) != 5) {
    x[4]
  }
})

# Some common structure of plots
comm <- ggplot(samples, aes(RNAseq, microbiota)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  theme(plot.title = element_text(hjust = 0.5))

samples <- cbind(samples, labels)
samples$Time <- factor(samples$Time, levels(as.factor(samples$Time))[c(1, 2, 4, 5, 3, 6, 7, 8)])
for (p in seq_along(levels(samples$Time))) {
  a <- comm +
    geom_text(aes(color = ID, label = ID)) +
    ggtitle(paste0("Samples by time")) +
    guides(col = guide_legend(title = "Patient")) +
    scale_color_manual(values = colors) +
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

for (p in seq_along(levels(samples$ID))) {
  a <- comm +
    geom_text(aes(color = ID, label = ifelse(!is.na(labels),
                                             paste(Time, labels, sep = "_"),
                                             as.character(Time)
    ))) +
    ggtitle(paste0("Samples by patient")) +
    guides(col = guide_legend(title = "Patient")) +
    scale_color_manual(values = colors) +
    facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = p)
  print(a)
}

comm +
  geom_text(aes(
    color = ID,
    label = ifelse(!is.na(labels),
                   paste(Time, labels, sep = "_"),
                   as.character(Time)
    )
  )) +
  guides(col = guide_legend(title = "Patient")) +
  scale_color_manual(values = colors)


comm +
  geom_text(aes(
    color = HSCT_responder,
    label = ifelse(!is.na(labels),
                   paste(Time, labels, sep = "_"),
                   as.character(Time)
    )
  )) +
  guides(col = guide_legend(title = "Responders"))

comm +
  geom_text(aes(
    color = Endoscopic_Activity,
    label = ifelse(!is.na(labels),
                   paste(ID, labels, sep = "_"),
                   as.character(ID)
    )
  )) +
  guides(col = guide_legend(title = "Endoscopic Activity"))

comm +
  geom_text(aes(color = Time, label = ifelse(!is.na(labels),
                                             paste(ID, labels, sep = "_"),
                                             as.character(ID)
  ))) +
  guides(col = guide_legend(title = "Time")) +
  scale_color_viridis_d()

comm +
  geom_text(aes(color = SESCD_local, label = ifelse(!is.na(labels),
                                             paste(ID, labels, sep = "_"),
                                             as.character(ID)
  ))) +
  guides(col = guide_legend(title = "SESCD (local)")) +
  scale_color_viridis_c()

variables <- data.frame(
  Origin = rep(names(A), sapply(A, ncol)),
  comp1 = unlist(sapply(
    sgcca.centroid$a,
    function(x) {
      x[, 1]
    }
  )),
  comp2 = unlist(sapply(
    sgcca.centroid$a,
    function(x) {
      x[, 2]
    }
  )),
  var = unlist(sapply(A, function(x){seq_len(ncol(x))}), use.names = FALSE)
)

# Remove the variables that in both components are 0
keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin == "16S"])
# keepComp1_metadata <- mean(abs(variables$comp1)[variables$Origin == "metadata"])

keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin == "16S"])
# keepComp2_metadata <- mean(abs(variables$comp2)[variables$Origin == "metadata"])

keepComp1 <- c(
  abs(variables$comp1[variables$Origin == "RNAseq"]) > keepComp1RNAseq,
  abs(variables$comp1[variables$Origin == "16S"]) > keepComp1_16S
  # abs(variables$comp1[variables$Origin == "metadata"]) > keepComp1_metadata
)
keepComp2 <- c(
  abs(variables$comp2[variables$Origin == "RNAseq"]) > keepComp2RNAseq,
  abs(variables$comp2[variables$Origin == "16S"]) > keepComp2_16S
  # abs(variables$comp2[variables$Origin == "metadata"]) > keepComp2_metadata
)

subVariables <- variables[keepComp1 & keepComp2, ]

ggplot(subVariables, aes(comp1, comp2), color = Origin) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
  geom_text(aes(color = Origin, label = var)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  coord_cartesian() +
  ggtitle(
    "Variables important for the first two components",
    subtitle = "Integrating stools and mucosa samples"
  )

rnaseq_i <- subVariables$var[subVariables$Origin == "RNAseq"]
if (length(rnaseq_i) >= 2) {
  pr <- prcomp(t(expr[rnaseq_i, ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
    ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
    ggtitle("RNAseq PCA from the important variables") +
    guides(col = guide_legend(title = "Responders"))
}

micro_i <- subVariables$var[subVariables$Origin == "16S"]
if (length(micro_i) >= 2) {
  pr <- prcomp(t(otus_table_i[micro_i, ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"] * 100)) +
    ylab(paste("PC2", prS$importance[2, "PC2"] * 100)) +
    ggtitle("16S PCA from the important variables") +
    guides(col = guide_legend(title = "Responders"))
}
# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x) {
  x[, 1]
})
variables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {
  x[, 2]
})
variables_weight(comp2)

l <- looIndex(size(A))
loo_model <- function(x, model){
  
  RGCCA::sgcca(A = subsetData(A, x),
               C = model, 
               scheme = "centroid", 
               verbose = FALSE, c1 = shrinkage
  )
}

result.out <- lapply(l, loo_model, model = model3)
saveRDS(result.out, "loo-model3.RDS")

# Bootstrap of sgcca
boot <- boot_sgcca(A, C, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model3.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()

