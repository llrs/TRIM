cd <- setwd("..")

rna <- "intestinal_RNAseq"

source("helper_functions.R")

# Read the intestinal RNAseq table

# Load the input data
expr <- read.delim(file.path(rna, "table.counts.results"), check.names = FALSE)

# Clean the RNAseq table
expr <- expr[rowSums(expr) != 0, ]

file_meta_r <- file.path(rna, "20171113_metadata.csv")
meta_r <- read.table(file_meta_r, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = ";", 
                     na.strings = c("NA", ""))
colnames(meta_r) <- meta_r[1, ]
meta_r <- meta_r[-1, ]

setwd(cd)


# Clean the metadata
meta_r <- meta_r[, apply(meta_r, 2, function(x){length(unique(x)) != 1})]
meta_r$ID <- meta_r$Patient_ID
meta_r$ID[meta_r$Patient_ID %in% c("15", "23")] <- "15/23"
meta_r$ID[meta_r$Patient_ID %in% c("33", "36")] <- "33/36"
meta_r$ID[meta_r$Patient_ID %in% c("29", "35")] <- "29/35"
meta_r$ID <- as.factor(meta_r$ID)


# Define colors
colors_i <- colors
names(colors_i) <- unique(meta_r$ID)

metadb <- meta_r

# Prepare the metadata
cols <- c(4, 7, 8, 13, 14, 15, 16, 17, 20, 21)
for (col in cols){
  metadb[, col] <- as.factor(metadb[, col])
}

for (col in cols){
  metadb[, col] <- as.numeric(metadb[, col])
}
metadb <- metadb[, cols]


metadb <- metadb[meta_r$`Sample Name_RNA` %in% colnames(expr), ]
metadb <- metadb[match(colnames(expr), meta_r$`Sample Name_RNA`), ]
metadb[is.na(metadb)] <- 0

meta_r <- meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ]
meta_r <- meta_r[match(colnames(expr), meta_r$`Sample Name_RNA`), ]

##### RGCCA #####
A <- list(intestinal = t(expr), meta = metadb)
C <- matrix(0, ncol = length(A), nrow = length(A), 
            dimnames = list(names(A), names(A)))
C <- subSymm(C, "intestinal", "meta", 1)

# Keep the covariance between them
shrinkage <- rep(1, length(A))

ncomp <- c(2, 2)

sgcca.centroid <-  sgcca(A, C, c1 = shrinkage,
                         ncomp = ncomp,
                         scheme = "centroid",
                         scale = TRUE,
                         verbose = FALSE)
names(sgcca.centroid$Y) <- names(A)
names(sgcca.centroid$a) <- names(A)
names(sgcca.centroid$astar) <- names(A)


PCA <- cbind(sgcca.centroid$Y[[1]], meta_r)
pdf(paste0("Figures/", today, "_plots.pdf"))

ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = Treatment)) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = HSCT_responder)) +
  guides(col = guide_legend(title="Responders")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = IBD)) +
  guides(col = guide_legend(title="Disease")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = ID)) +
  guides(col = guide_legend(title="Patient")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = CD_Aftected_area))  +
  guides(col = guide_legend(title="Afected area")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = Exact_location)) +
  guides(col = guide_legend(title="Exact location")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = Involved_Healthy)) +
  guides(col = guide_legend(title="Healthy")) +
  ggtitle("Intestinal RNAseq PCA-like")
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = Time)) +
  ggtitle("Intestinal RNAseq PCA-like")
dev.off()
