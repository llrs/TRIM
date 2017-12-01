cd <- setwd("..")

intestinal <- "intestinal_16S"

source("helper_functions.R")

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)

# Clean the otu table
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]



file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
setwd(cd)


# Clean the metadata
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i <- meta_i[, apply(meta_i, 2, function(x){length(unique(x)) != 1})]
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i$ID <- meta_i$Patient_ID
meta_i$ID[meta_i$Patient_ID %in% c("15", "23")] <- "15/23"
meta_i$ID[meta_i$Patient_ID %in% c("33", "36")] <- "33/36"
meta_i$ID[meta_i$Patient_ID %in% c("29", "35")] <- "29/35"
meta_i$ID <- as.factor(meta_i$ID)


# Define colors
colors_i <- colors
names(colors_i) <- unique(meta_i$ID)

metadb <- meta_i
# Prepare the metadata
for (col in seq_len(ncol(metadb[, c(3:11, 13)]))){
  metadb[, c(3:11, 13)][, col] <- as.factor(metadb[, c(3:11, 13)][, col])
}

for (col in seq_len(ncol(metadb[, c(3:10, 13)]))){
  levels(metadb[, c(3:10, 13)][, col]) <- seq_along(levels(metadb[, c(3:10, 13)][, col]))
}
metadb <- apply(metadb[, c(3:11, 13)], 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Select variables that explain the disease
keepCol <- c("Time", "Age", "IBD", "Active_area", "CD_Aftected_area", 
             "Involved_Healthy", "Endoscopic_Activity", "ID")

##### RGCCA #####
A <- list(intestinal = t(otus_table_i), meta = metadb[, keepCol])
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


PCA <- cbind(sgcca.centroid$Y[[1]], meta_i)
pdf(paste0("Figures/", today, "_plots.pdf"))
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = Treatment))
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = HSCT_responder))
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = IBD))
ggplot(PCA) + 
  geom_point(aes(comp1, comp2, col = ID))
dev.off()