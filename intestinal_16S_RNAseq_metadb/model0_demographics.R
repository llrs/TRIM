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
model0l <- subSymm(model, "16S", "Demographics", 1)
model0l <- subSymm(model0l, "RNAseq", "Demographics", 1)
model0li <- subSymm(model0l, 1, 1, 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.25670333, 0, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

sgcca.centroid <- sgcca(
  A, C = model0l, c1 = shrinkage,
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
  "microbiota" = sgcca.centroid$Y[["16S"]][, 1],
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
  geom_text(aes(
    color = Exact_location,
    label = ifelse(!is.na(labels),
                   paste(ID, labels, sep = "_"),
                   as.character(ID)
    )
  )) +
  guides(col = guide_legend(title = "Location"))

comm +
  geom_text(aes(color = Time, label = ifelse(!is.na(labels),
                                             paste(ID, labels, sep = "_"),
                                             as.character(ID)
  ))) +
  guides(col = guide_legend(title = "Time"))

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
  ))
)
variables$var <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
variables$var <- gsub("^RNAseq\\.(ENSG.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(ENSG.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))
variables$var <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))

# Remove the variables that in both components are 0
keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin == "16S"])
keepComp1_metadata <- mean(abs(variables$comp1)[variables$Origin == "metadata"])

keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin == "16S"])
keepComp2_metadata <- mean(abs(variables$comp2)[variables$Origin == "metadata"])

keepComp1 <- c(
  abs(variables$comp1[variables$Origin == "RNAseq"]) > keepComp1RNAseq,
  abs(variables$comp1[variables$Origin == "16S"]) > keepComp1_16S,
  abs(variables$comp1[variables$Origin == "metadata"]) > keepComp1_metadata
)
keepComp2 <- c(
  abs(variables$comp2[variables$Origin == "RNAseq"]) > keepComp2RNAseq,
  abs(variables$comp2[variables$Origin == "16S"]) > keepComp2_16S,
  abs(variables$comp2[variables$Origin == "metadata"]) > keepComp2_metadata
)

subVariables <- variables[keepComp1 & keepComp2, ]

ggplot(subVariables, aes(comp1, comp2), color = Origin) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
  geom_text(aes(color = Origin, label = trimVer(var))) +
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


# Validate ####

l <- looIndex(size(A))
result.out <- lapply(l, function(x){
  
  RGCCA::sgcca(A = subsetData(A, x),
               C = model0l, 
               scheme = "centroid", 
               verbose = FALSE, c1 = shrinkage
  )}) # SGCCA of the selected model leaving one sample each time out of order.
saveRDS(result.out, "loo-model0l.RDS")
result.out <- lapply(l, function(x){
  
  RGCCA::sgcca(A = subsetData(A, x),
               C = model0li, 
               scheme = "centroid", 
               verbose = FALSE, c1 = shrinkage
  )}) # SGCCA of the selected model leaving one sample each time out of order.
saveRDS(result.out, "loo-model0li.RDS")


# Bootstrap of sgcca
boot <- boot_sgcca(A, model0l, shrinkage, 1000)

saveRDS(boot, file = "bootstrap.RDS")
boot <- boot_sgcca(A, model0li, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model0l.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()
