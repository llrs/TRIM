cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

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

# Correct metadata
meta_i <- meta_i_norm(meta_i)
meta_r <- meta_r_norm(meta_r)

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Filter expression
expr <- norm_RNAseq(expr)

#' Remove otus not present in that set of samples
otus_table_i <- otus_table_i[rowSums(otus_table_i) > 0, ]

# Prepare input for the sgcca function
A <- list("RNAseq" = t(expr), "16S" = t(otus_table_i))
A <- sapply(A, function(x){
  x[, apply(x, 2, sd) != 0]
}, simplify = FALSE)
save(A, meta_i, file = "TRIM.RData")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "RNAseq", 1)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.1, 0) # We guess a 0.5
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the threshold  allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))

ncomp <- c(2, 2)

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

sgcca.factorial <- sgcca(
  A, C, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "factorial",
  scale = TRUE,
  verbose = FALSE
)
names(sgcca.factorial$Y) <- names(A)
names(sgcca.factorial$a) <- names(A)
names(sgcca.factorial$astar) <- names(A)
names(sgcca.factorial$AVE$AVE_X) <- names(A)
sgcca.factorial$AVE$AVE_X <- simplify2array(sgcca.factorial$AVE$AVE_X)

sgcca.horst <- sgcca(
  A, C, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "horst",
  scale = TRUE,
  verbose = FALSE
)
names(sgcca.horst$Y) <- names(A)
names(sgcca.horst$a) <- names(A)
names(sgcca.horst$astar) <- names(A)
names(sgcca.horst$AVE$AVE_X) <- names(A)
sgcca.horst$AVE$AVE_X <- simplify2array(sgcca.factorial$AVE$AVE_X)

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)


samples <- data.frame(
  Stools = sgcca.centroid$Y[[1]][, 1],
  Intestinal = sgcca.centroid$Y[[2]][, 1]
)
lmt <- lm(Intestinal ~ Stools, data = samples)

if (lmt$coefficients[2] > 0) {
  d <- c(1, 1)
} else if (lmt$coefficients[2] < 0) {
  d <- c(-1, -1)
}


dist <- apply(samples, 1, dist2d, d = d)
# Colors for the plots
names(colors) <- unique(meta_r$ID)

samples <- cbind(samples, meta_r, "dist" = dist)
samples$Patient_ID <- as.factor(samples$Patient_ID)
samples$Sample_Code <- as.character(samples$Sample_Code)

pdf(paste0("Figures/", today, "_RGCCA_plots.pdf"))

# Plot if the coherence between samples has a specific pattern
ggplot(samples) +
  geom_point(aes(Patient_ID, log10(dist), col = Involved_Healthy)) +
  facet_grid(~ Time)

hist(samples$dist)

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

samples <- cbind(samples, labels)
samples$Time <- factor(samples$Time, levels(as.factor(samples$Time))[c(1, 2, 4, 5, 3, 6, 7, 8)])
for (p in seq_along(levels(samples$Time))) {
  a <- ggplot(samples, aes(Stools, Intestinal)) +
    geom_text(aes(color = ID, label = ID)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Patient")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    geom_abline(intercept = 0, slope = d[1], linetype = 2) +
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

for (p in seq_along(levels(samples$ID))) {
  a <- ggplot(samples, aes(Stools, Intestinal)) +
    geom_text(aes(color = ID, label = ifelse(!is.na(labels),
      paste(Time, labels, sep = "_"),
      as.character(Time)
    ))) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by patient")) +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Patient")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    geom_abline(intercept = 0, slope = d[1], linetype = 2) +
    facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = p)
  print(a)
}
ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(
    color = ID,
    label = ifelse(!is.na(labels),
      paste(Time, labels, sep = "_"),
      as.character(Time)
    )
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Patient")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors) +
  geom_abline(intercept = 0, slope = d[1], linetype = 2)


ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(
    color = HSCT_responder,
    label = ifelse(!is.na(labels),
      paste(ID, labels, sep = "_"),
      as.character(Time)
    )
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Responders")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_abline(intercept = 0, slope = d[1], linetype = 2)


ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(
    color = Endoscopic_Activity,
    label = ifelse(!is.na(labels),
      paste(ID, labels, sep = "_"),
      as.character(ID)
    )
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Endoscopic Activity")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_abline(intercept = 0, slope = d[1], linetype = 2)

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color = Time, label = ifelse(!is.na(labels),
    paste(ID, labels, sep = "_"),
    as.character(ID)
  ))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Time")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_abline(intercept = 0, slope = d[1], linetype = 2)

variables <- data.frame(
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
  Origin = rep(names(A), sapply(A, ncol))
)
variables$var <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
variables$var <- gsub("^RNAseq\\.(ENSG.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(ENSG.*)$", "\\1", rownames(variables))

# Remove the variables that in both components are 0
keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin != "RNAseq"])

keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin != "RNAseq"])

keepComp1 <- c(
  variables$comp1[variables$Origin == "RNAseq"] > keepComp1RNAseq,
  variables$comp1[variables$Origin != "RNAseq"] > keepComp1_16S
)
keepComp2 <- c(
  variables$comp2[variables$Origin == "RNAseq"] > keepComp2RNAseq,
  variables$comp2[variables$Origin != "RNAseq"] > keepComp2_16S
)
subVariables <- variables[keepComp1 | keepComp2, ]

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

save(sgcca.centroid, file = "sgcca.RData")

# To calculate the conficence interval on selecting the variable
# this interval should reduce as we fit a better model/relationship
# Bootstrap of sgcca
STAB <- boot_sgcca(A, C, shrinkage, 1000)

save(STAB, file = "bootstrap.RData")

# Evaluate the boostrap effect and plot
boot_evaluate(STAB)

dev.off()
