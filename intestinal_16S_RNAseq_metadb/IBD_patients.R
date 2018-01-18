cd <- setwd("..")
library("fgsea")

# Load the helper file
source("helper_functions.R")
source("helper_RGCCA.R")

intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"


# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))


# Load the input data
expr <- read.delim(file.path(rna, "table.counts.results"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
file_meta_r <- file.path(rna, "111217_metadata.csv")
meta_r <- read.table(file_meta_r, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = ";",
                     na.strings = c(NA, ""), header = TRUE, dec = c(",", "."))

setwd(cd)


# Filter only the IBD patients
meta_r <- meta_r[meta_r$HSCT_responder != "C", ]

# Find the samples that we have microbiota and expression
int <- intersect(meta_r$Sample_Code_uDNA[!is.na(meta_r$Sample_Code_uDNA) & 
                                           !is.na(meta_r$`Sample Name_RNA`)],
                 meta_i$Sample_Code)
# table(sapply(strsplit(int, "_"), "[", 1))

# Correct the metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

meta_i <- meta_i[meta_i$Sample_Code %in% int, ]
meta_r <- meta_r[meta_r$Sample_Code_uDNA %in% int, ]
meta_r <- meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Match the labels and order to append the id
meta_i <- meta_i[match(meta_r$Sample_Code_uDNA, meta_i$Sample_Code), ]
meta_r$`Sample Name_Code` <- gsub("([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1", rownames(meta_i))

colnames(otus_table_i) <- gsub("([0-9]{2,3}\\.B[0-9]+)\\..+", "\\1", 
                               colnames(otus_table_i))

# Subset expression and outs
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$`Sample Name_Code`]

# Subset if all the rows are 0 and if sd is 0
otus_table_i <- otus_table_i[apply(otus_table_i, 1, sd) != 0, ] 

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them 
metadb <- meta_r
keepCol <- sapply(metadb, is.factor)
nam <- c("Involved_Healthy", 
         "Active_area", "AGE_SAMPLE", "Transplant", "ID", 
         "Exact_location", "Endoscopic_Activity", "Treatment", "SEX")
keepCol <- keepCol[nam]
keepCol[nam] <- TRUE
for (col in names(keepCol)){
  if (class(metadb[, col]) == "character") {
    metadb[, col] <- as.factor(metadb[, col])
    levels(metadb[, col]) <- seq_along(levels(metadb[, col]))
  } else if (class(metadb[, col]) == "factor") {
    levels(metadb[, col]) <- seq_along(levels(metadb[, col]))
  } else if (class(metadb[, col]) == "numeric") {
    next
  }
}
metadb <- metadb[, names(keepCol)]

# Set metadb with a sigle variable with several options
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Filter expression
expr <- norm_RNAseq(expr)

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "metadata" = metadb)
saveRDS(A, file = "TRIM_IBD.RDS")

# The design
C <- matrix(0, ncol = length(A), nrow = length(A), 
            dimnames = list(names(A), names(A)))
C <- subSymm(C, "16S", "metadata", 1)
C <- subSymm(C, "RNAseq", "metadata", 1)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.122747, 0, 1) # We guess a 0.5 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x){1/sqrt(ncol(x))}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- c(2, 2, 2)

sgcca.centroid <-  sgcca(A, C, c1 = shrinkage,
                         ncomp = ncomp,
                         scheme = "centroid",
                         scale = TRUE,
                         verbose = FALSE)
names(sgcca.centroid$Y) <- names(A)
names(sgcca.centroid$a) <- names(A)
names(sgcca.centroid$astar) <- names(A)

sgcca.factorial <-  sgcca(A, C, c1 = shrinkage,
                          ncomp = ncomp,
                          scheme = "factorial",
                          scale = TRUE,
                          verbose = FALSE)
names(sgcca.factorial$Y) <- names(A)
names(sgcca.factorial$a) <- names(A)
names(sgcca.factorial$astar) <- names(A)

sgcca.horst <-  sgcca(A, C, c1 = shrinkage,
                      ncomp = ncomp,
                      scheme = "horst",
                      scale = TRUE,
                      verbose = FALSE)
names(sgcca.horst$Y) <- names(A)
names(sgcca.horst$a) <- names(A)
names(sgcca.horst$astar) <- names(A)

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst, 
# sgcca.factorial = sgcca.factorial)
save(sgcca.centroid, file = "IBD.RData")

samples <- data.frame("RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
                      "microbiota" = sgcca.centroid$Y[["16S"]][, 1])


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

km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
plot(samples[, c("RNAseq", "microbiota")], col = km$cluster)

## Plotting results ####
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, meta_r)
samples$Patient_ID <- as.factor(samples$Patient_ID)
samples$Sample_Code <- as.character(samples$Sample_Code)

pdf(paste0("Figures/", today, "_RGCCA_plots_IBD.pdf"))

# Labels of the samples
label <- strsplit(as.character(samples$`Sample Name_RNA`), split = "-")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    x[5]
  }
  else if (length(x) != 5) {
    x[4]
  }
})

samples <- cbind(samples, labels)
samples$Time <- factor(samples$Time, levels(as.factor(samples$Time))[c(1, 2, 4, 5, 3, 6, 7, 8)])
for (p in seq_along(levels(samples$Time))){
  a <- ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color =  ID, label = ID)) + 
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) + 
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title="Patient")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) + 
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

for (p in seq_along(levels(samples$ID))){
  a <- ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color =  ID, label = ifelse(!is.na(labels), 
                                              paste(Time, labels, sep = "_"),
                                              as.character(Time)))) + 
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by patient")) + 
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title="Patient")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) + 
    facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = p)
  print(a)
}
ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(color =  ID, 
                label = ifelse(!is.na(labels), 
                               paste(Time, labels, sep = "_"),
                               as.character(Time)))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  # xlab("RNAseq (component 1)") +
  # ylab("16S (component 1)") +
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors)


ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(color = HSCT_responder, 
                label = paste(Time, labels, sep = "_"))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title="Responders")) + 
  theme(plot.title = element_text(hjust = 0.5))


ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(color = Endoscopic_Activity , 
                label = ifelse(!is.na(labels), 
                               paste(ID, labels, sep = "_"),
                               as.character(ID)))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Endoscopic Activity")) + 
  theme(plot.title = element_text(hjust = 0.5)) 

ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(color = Time , label = ifelse(!is.na(labels), 
                                              paste(ID, labels, sep = "_"),
                                              as.character(ID)))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Time")) + 
  theme(plot.title = element_text(hjust = 0.5))

variables <- data.frame(Origin = rep(names(A), sapply(A, ncol)),
                        comp1 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 1]})),
                        comp2 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 2]})))
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

keepComp1 <- c(abs(variables$comp1[variables$Origin == "RNAseq"]) > keepComp1RNAseq,
               abs(variables$comp1[variables$Origin == "16S"]) > keepComp1_16S,
               abs(variables$comp1[variables$Origin == "metadata"]) > keepComp1_metadata)
keepComp2 <- c(abs(variables$comp2[variables$Origin == "RNAseq"]) > keepComp2RNAseq,
               abs(variables$comp2[variables$Origin == "16S"]) > keepComp2_16S,
               abs(variables$comp2[variables$Origin == "metadata"]) > keepComp2_metadata)

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
  ggtitle("Variables important for the first two components", 
          subtitle = "Integrating stools and mucosa samples")

rnaseq_i <- subVariables$var[subVariables$Origin == "RNAseq"]
if (length(rnaseq_i) >= 2) {
  pr <- prcomp(t(expr[rnaseq_i, ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"]*100)) +
    ylab(paste("PC2", prS$importance[2, "PC2"]*100)) +
    ggtitle("RNAseq PCA from the important variables")
}

micro_i <- subVariables$var[subVariables$Origin == "16S"]
if (length(micro_i) >= 2) {
  pr <- prcomp(t(otus_table_i[micro_i, ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta_r$HSCT_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"]*100)) +
    ylab(paste("PC2", prS$importance[2, "PC2"]*100)) +
    ggtitle("16S PCA from the important variables")
}
# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x){x[, 1]})
Loadings <- unlist(comp1)
comp1 <- as.data.frame(Loadings)
comp1$Origin <- factor(rep(names(sgcca.centroid$a), 
                           lengths(sgcca.centroid$a)/2), 
                       levels = names(sgcca.centroid$a))
rownames(comp1) <- seq_len(nrow(comp1))
ggplot(comp1) +
  stat_density(aes(x = Loadings, 
                   y = ..scaled..,
                   fill = Origin), alpha = 0.5) +
  ggtitle("Importance of the otus of each data set") +
  ylab("Scaled density") +
  xlab("weight") +
  facet_grid(~Origin) + 
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x){x[, 2]})
Loadings <- unlist(comp2)
comp2 <- as.data.frame(Loadings)
comp2$Origin <- comp1$Origin
rownames(comp2) <- seq_len(nrow(comp2))
ggplot(comp2) +
  stat_density(aes(x = Loadings, y = ..scaled.., fill = Origin), alpha = 0.5) +
  ggtitle("Importance of each block variable", 
          subtitle = "Second component") +
  ylab("Scaled density") +
  xlab("weight") +
  facet_grid(~Origin) + 
  guides(fill = FALSE) 

# Bootstrap of sgcca 
STAB <- boot_sgcca(A, C, shrinkage, 1000)

save(STAB, file = "bootstrap_IBD.RData")

# Evaluate the boostrap effect and plot 
boot_evaluate(STAB)


dev.off()