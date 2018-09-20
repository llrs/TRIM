cd <- setwd("..")
library("ggforce")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

library("fgsea")

intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Load the input data
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(cd)

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- colnames(expr)[rev(position)]
colnames(expr) <- toupper(colnames(expr))
#To match metadata
colnames(expr) <- gsub("16-TM29", "16-TM30", colnames(expr)) 

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes=TRUE, log = TRUE)

# Filter expression
expr <- norm_RNAseq(expr_norm)

# Normalize OTUS
library("metagenomeSeq")
MR_i <- newMRexperiment(
  otus_table_i, 
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i[rownames(otus_table_i), ]))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus_table_i <- MRcounts(MR_i, norm = TRUE, log = TRUE)

# Subset if all the rows are 0 and if sd is 0
otus_table_i <- otus_table_i[apply(otus_table_i, 1, sd) != 0, ]

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

Localization <- model_RGCCA(metadb, c("Exact_location")) # With SESCD local it increase the AVE_inner
Time <- model_RGCCA(metadb, c("AgeDiag", "AGE_SAMPLE", "Transplant"))
Demographics <- model_RGCCA(metadb, c("ID","SEX", "Surgery", "Treatment"))
Time[is.na(Time)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "Demographics" = Demographics,
          "Location" = Localization, "Time" = Time)

stopifnot(length(unique(vapply(A, nrow, numeric(1L)))) == 1)

# The design of model 3
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "Demographics", 1)
C <- subSymm(C, "16S", "Time", 1)
C <- subSymm(C, "16S", "Location", 1)
C <- subSymm(C, "RNAseq", "Demographics", 1)
C <- subSymm(C, "RNAseq", "Time", 1)
C <- subSymm(C, "RNAseq", "Location", 1)
C <- subSymm(C, "RNAseq", "16S", 1)
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
  A, C, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
names(sgcca.centroid$Y) <- names(A)
names(sgcca.centroid$a) <- names(A)
names(sgcca.centroid$astar) <- names(A)
names(sgcca.centroid$AVE$AVE_X) <- names(A)
sgcca.centroid$AVE$AVE_X <- simplify2array(sgcca.centroid$AVE$AVE_X)
sgcca.centroid$AVE

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)
saveRDS(sgcca.centroid, file = "sgcca_model3.RDS")

designs <- weight_design(3, 5)
designs_rank <- sapply(designs, function(x){
  rankMatrix(x)
})

design_boot <- sapply(designs[designs_rank == 5], function(x) {
  sgcca.centroid <- RGCCA::sgcca(
    A, x, c1 = shrinkage,
    ncomp = ncomp,
    scheme = "centroid",
    scale = TRUE,
    verbose = FALSE
  )
  names(sgcca.centroid$Y) <- names(A)
  names(sgcca.centroid$a) <- names(A)
  names(sgcca.centroid$astar) <- names(A)
  names(sgcca.centroid$AVE$AVE_X) <- names(A)
  sgcca.centroid$AVE$AVE_X <- simplify2array(sgcca.centroid$AVE$AVE_X)
  sgcca.centroid$AVE[c("AVE_inner", "AVE_outer")]
})
saveRDS(design_boot, "designs_boot_model3.RDS")

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
  guides(col = guide_legend(title = "Time"))

comm +
  geom_text(aes(color = SESCD_local, label = ifelse(!is.na(labels),
                                             paste(ID, labels, sep = "_"),
                                             as.character(ID)
  ))) +
  guides(col = guide_legend(title = "SESCD (local)"))

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

# Bootstrap of sgcca
boot <- boot_sgcca(A, C, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model3.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()
