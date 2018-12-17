cd <- setwd("..")

stool <- "stools_16S"
today <- format(Sys.time(), "%Y%m%d")
library("integration")


# Read the stools OTUs
otus_table_s <- read.delim(
  file.path(stool, "OTUs-Table-refined-stools.tab"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(
  file_meta_s, check.names = FALSE, row.names = 1,
  stringsAsFactors = FALSE
)
setwd(cd)


# Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x) {
  length(unique(x)) != 1
})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)

# Define colors
colors_s <- colors
names(colors_s) <- unique(meta_s$ID)

# Select variables that explain the disease
keepCol <- c("Time", "Age", "ID")
metadb <- meta_s[, keepCol]

# Prepare the metadata
for (col in seq_len(ncol(metadb))) {
  metadb[, col] <- as.factor(metadb[, col])
  levels(metadb[, col]) <- seq_along(levels(metadb[, col]))
}
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

##### RGCCA #####
A <- list(stool = t(otus_table_s), meta = metadb)
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "stool", "meta", 1)

# Keep the covariance between them
shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

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


PCA <- cbind(sgcca.centroid$Y[[1]], meta_s)
pdf(paste0("Figures/", today, "_plots.pdf"))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = Treatment))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = ID))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = Endoscopic_Activity))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = HSCT_responder))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = Age))

dev.off()
