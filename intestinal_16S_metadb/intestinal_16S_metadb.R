cd <- setwd("..")

intestinal <- "intestinal_16S"

today <- format(Sys.time(), "%Y%m%d")
library("integration")

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)

# Clean the otu table
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]



file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
setwd(cd)


# Clean the metadata
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i <- meta_i[, apply(meta_i, 2, function(x) {
  length(unique(x)) != 1
})]
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i$ID <- meta_i$Patient_ID
meta_i$ID[meta_i$Patient_ID %in% c("15", "23")] <- "15/23"
meta_i$ID[meta_i$Patient_ID %in% c("33", "36")] <- "33/36"
meta_i$ID[meta_i$Patient_ID %in% c("29", "35")] <- "29/35"
meta_i$ID <- as.factor(meta_i$ID)

# Pre transplant
meta_i$Transplant <- "Post" #
meta_i$Transplant[meta_i$Patient_ID %in% c("15", "33", "29")] <- "Pre"
meta_i$Transplant[grep("^C", meta_i$Patient_ID)] <- NA

# There is a mislabeling on those tubes, we don't know which is which
meta_i$CD_Aftected_area[meta_i$Sample_Code == "22_T52_T_DM_III"] <- NA


# Define colors
colors_i <- colors
names(colors_i) <- unique(meta_i$ID)

metadb <- meta_i
cols <- c(3:10, 13, 14)
# Prepare the metadata
for (col in cols) {
  metadb[, col] <- as.factor(metadb[, col])
  levels(metadb[, col]) <- seq_along(levels(metadb[, col]))
}

cols <- c(cols, 11)
metadb <- apply(metadb[, cols], 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

##### RGCCA #####
A <- list(intestinal = t(otus_table_i), meta = metadb)
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "intestinal", "meta", 1)

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
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))

PCA <- cbind(sgcca.centroid$Y[[1]], meta_i)
pdf(paste0("Figures/", today, "_plots.pdf"))
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = Treatment)) +
  ggtitle("Intestinal 16S PCA-like")
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = HSCT_responder)) +
  guides(col = guide_legend(title = "Responders")) +
  ggtitle("Intestinal 16S PCA-like")
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = IBD)) +
  guides(col = guide_legend(title = "Disease")) +
  ggtitle("Intestinal 16S PCA-like")
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = ID)) +
  guides(col = guide_legend(title = "Patient")) +
  ggtitle("Intestinal 16S PCA-like")
ggplot(PCA) +
  geom_point(aes(comp1, comp2, col = Time)) +
  ggtitle("Intestinal 16S PCA-like")
dev.off()
