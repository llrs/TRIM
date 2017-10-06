intestinal <- "intestinal_16S"
stool <- "stools_16S"
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

#' Clean and prepare the data from IMNGS
#' 
#' Divides the taxonomy into a new matrix for each otu
#' 
#' @param taxonomy Last column of files a string ; separated with domain, 
#' phylum, vlass, order, family, genus and species.
#' @param otus The name of the rows
#' 
#' @return 
#' A matrix with the taxonomic information ready for the package phylo
taxonomy <- function(taxonomy, otus){
  taxonomy <- sapply(taxonomy, strsplit, split = ";")
  names(taxonomy) <- otus
  otus_tax <- t(sapply(taxonomy, '[', seq(max(sapply(taxonomy, length)))))
  colnames(otus_tax) <- c("Domain", "Phylum", "Class", "Order", 
                          "Family", "Genus", "Species")
  # Remove spaces
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "\\s", replacement = "")
  otus_tax[otus_tax == ""] <- NA # Remove empty cells
  otus_tax
}

otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Read for the stools
otus_table_s <- read.delim(file.path(stool, "OTUs-Table-refined-stools.tab"), 
                           stringsAsFactors = FALSE, row.names = 1,
                           check.names = FALSE)
tax_s <- otus_table_s[, ncol(otus_table_s)]
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]
otus_tax_s <- taxonomy(tax_s, rownames(otus_table_s))

# Check the taxonomy
# https://stackoverflow.com/q/7943695/2886003
fastercheck <- function(x,matrix){
  nc <- ncol(matrix)
  rec.check <- function(r,i,id){
    id[id] <- matrix[id,i] %in% r[i]
    if(i<nc & any(id)) rec.check(r,i+1,id) else any(id)
  }
  apply(x,1,rec.check,1,rep(TRUE,nrow(matrix)))
}
eq <- sum(fastercheck(as.matrix(otus_tax_i), as.matrix(otus_tax_s)))
(dice <- 2*eq/(nrow(otus_tax_i)+nrow(otus_tax_s)))

# Read the metadata for each type of sample
files_s <- list.files(path = stool, pattern = ".txt", full.names = TRUE)
meta_s <- read.delim(files_s[1], check.names = FALSE, row.names = 1, 
                  stringsAsFactors = FALSE)
files_i <- list.files(path = intestinal, pattern = ".txt", full.names = TRUE)
meta_i <- read.delim(files_i[1], row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)

# Reorder by patient and time
meta_s <- meta_s[order(meta_s$Patient_ID, meta_s$Time), ]
meta_i <- meta_i[order(meta_i$Patient_ID, meta_i$Time), ]

# Find common patients
comPatient <- intersect(meta_i$Patient_ID, meta_s$Patient_ID)
comTime <- intersect(meta_i$Time, meta_s$Time)

# Keep only the common patients and times
com_meta_i <- meta_i[meta_i$Patient_ID %in% comPatient & meta_i$Time %in% comTime, ]
com_meta_s <- meta_s[meta_s$Patient_ID %in% comPatient & meta_s$Time %in% comTime, ]

# Samples per patient and Time
tab_s <- table(com_meta_s$Patient_ID, com_meta_s$Time)
tab_i <- table(com_meta_i$Patient_ID, com_meta_i$Time)

# Compare if we have more intestinal samples than stools for the same time and
# patient
moreS <- tab_i > tab_s & tab_i != 0 & tab_s != 0
QmoreS <- tab_i - tab_s

# Delete rows of patients which are not in common between the datasets
keep_i <- rownames(meta_i)[meta_i$Patient_ID %in% comPatient & 
                             meta_i$Time %in% comTime]
keep_s <- rownames(meta_s)[meta_s$Patient_ID %in% comPatient & 
                             meta_s$Time %in% comTime]
com_otus_table_i <- otus_table_i[, keep_i]
com_otus_table_s <- otus_table_s[, keep_s]

# Extract the ids of the data we need to add or remove
library("reshape2")
df <- melt(QmoreS, varnames = c("Patient_ID", "Time"))
df <- df[df$value != 0, ]

# Remove the stools that should be removed because we don't have intestinal 
# samples
remove_s <- df[df$value <0, ]
rownames(remove_s) <- seq_len(nrow(remove_s))

remove_samples_s <- apply(remove_s, 1, function(x){
  rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]]
})

com_otus_table_s <- com_otus_table_s[, 
                                     !colnames(com_otus_table_s) %in% remove_samples_s]

# Add stool samples as much as need
add_s <- df[df$value >0, ]
rownames(add_s) <- seq_len(nrow(add_s))

remove_samples_s <- apply(add_s, 1, function(x){
  rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]]
})

com_otus_table_s <- com_otus_table_s[, 
                                     !colnames(com_otus_table_s) %in% remove_samples_s]


# # Integrate them
# library("RGCCA")
# A <- list(stools = otus_table_s, intestinal = otus_table_i, metadata = )
# C <- matrix(ncol = 2, nrow = 2)
# sgcca.glioma = sgcca(A, C, c1 = c(.071,.2, 1),
#                      ncomp = c(1, 1, 1),
#                      scheme = "centroid",
#                      scale = TRUE,
#                      verbose = FALSE)