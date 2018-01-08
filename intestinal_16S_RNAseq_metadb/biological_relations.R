# Test if the loadings of the genes are from a relevant pathways. 
wd <- setwd("..")

# Load the helper file
source("helper_functions.R")

intestinal <- "intestinal_16S"
stool <- "stools_16S"
rna <- "intestinal_RNAseq"
source("helper_functions.R")

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))


setwd(wd)

library("fgsea") # To test the distribution of the pathways
library("reactome.db") # To get the pathways
library("org.Hs.eg.db") # To get the names to entrez from the ensembl
library("data.table")

# Load previous data
load("sgcca.RData")
load("bootstrap.RData")

# Select the significant relationships ####
#' Two sided test
#' 
#' Test in a vector from a permutation if there is a relationship or not. 
#' Assumes that the distribution is symmetric around 0.
#' @param z Vector of of the permutations
#' @param y Value of the test
#' @return The p-value
two.sided <- function(y, z) {
  stopifnot(length(y) == 1)
  stopifnot(length(z) > 2)
  sum(abs(z) > abs(y), na.rm = TRUE)/length(z)
}

b <- STAB[["RNAseq"]]
d <- sgcca.centroid$a[["RNAseq"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!duplicated(b[, 1:10]), ]

pvalue <- numeric(ncol(b))
for (col in seq_len(ncol(b))) {
  pvalue[col] <- two.sided(d[col], b[, col])
}
names(pvalue) <- colnames(b)
fdr <- p.adjust(pvalue, "BH")

# Select those genes that are significant
significant <- names(fdr)[fdr < 0.05]


loadings <- sgcca.centroid$a[["RNAseq"]]

# Convert the ids to the entrezIds
ensemblID <- rownames(loadings)
ensemblID <- gsub("(.*)\\..*", "\\1", ensemblID)
rownames(loadings) <- gsub("(.*)\\..*", "\\1", rownames(loadings))
entrezID <- mapIds(org.Hs.eg.db, keys = ensemblID, keytype = "ENSEMBL", 
                   column = "ENTREZID")
comp1 <- loadings[, 1]
names(comp1) <- entrezID

# Extract the information of the pathways
genes2Pathways <- as.list(reactome.db::reactomeEXTID2PATHID)
pathways <- unlist(genes2Pathways, use.names = FALSE)
genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
paths2genes <- split(genes, pathways) # List of genes and the gene sets

# Subset to only human pathways
paths2genes <- paths2genes[grep("R-HSA-", names(paths2genes))]
path2genes <- c(path2genes, "significant" = significant)

# Compute the GSEA for the size effect
gseaSizeEffect <- fgsea(paths2genes, comp1, nperm = length(comp1))

# Get the name of the pathway
namesPaths <- select(reactome.db, keys = gseaSizeEffect$pathway, 
                     keytype = "PATHID", columns = "PATHNAME")
# Remove the homo sapiens part
namesPaths$PATHNAME <- gsub("Homo sapiens: (.*)", "\\1", namesPaths$PATHNAME)
# Add a column
gseaSizeEffect[ , namesPaths := namesPaths$PATHNAME]
# Order the dataframe by size effect 
data.table::setorder(gseaSizeEffect, -NES, padj, -size)
# Store the output
fwrite(gseaSizeEffect[padj < 0.05, ], file = "relevant_Effect_pathways.csv")

hist(fdr, main = "Histogram of FDR for RNAseq", xlab = "FDR")
abline(v = 0.05, col = "red")

# Compute the GSEA
gseaFDR <- fgsea(paths2genes, fdr, nperm = 5000)
# Get the name of the pathway
namesPaths <- select(reactome.db, keys = gseaFDR$pathway, 
                     keytype = "PATHID", columns = "PATHNAME")
# Remove the homo sapiens part
namesPaths$PATHNAME <- gsub("Homo sapiens: (.*)", "\\1", namesPaths$PATHNAME)
# Add a column
gseaFDR[ , namesPaths := namesPaths$PATHNAME]
# Order the dataframe by size effect 
data.table::setorder(gseaFDR, -NES, padj, -size)
# Store the output
fwrite(gseaFDR[padj < 0.05, ], file = "relevant_Related_pathways.csv")

sign_size <- intersect(gseaSizeEffect$pathway[gseaSizeEffect$padj < 0.05],
                       gseaFDR$pathway[gseaFDR$padj < 0.05])


##
b <- STAB[["16S"]]
d <- sgcca.centroid$a[["16S"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!duplicated(b[, 1:10]), ]

pvalue <- numeric(ncol(b))
for (col in seq_len(ncol(b))) {
  pvalue[col] <- m(d[col], b[, col])
}
names(pvalue) <- names(d)
fdr <- p.adjust(pvalue, "BH")
hist(fdr, main = "Histogram of FDR for RNAseq", xlab = "FDR")
abline(v = 0.05, col = "red")
sign_otus <- data.frame(loading = as.numeric(d[fdr < 0.05]), 
                        otus_tax_i[names(d[fdr < 0.05]), c("Genus", "Species")])
data.table::setDT(sign_otus)
fwrite(sign_otus, "relevant_otus.csv")
sign_otus[, .(mean = mean(loading), sd = sd(loading)), keyby = Genus]

## Names
b <- STAB[["metadata"]]
d <- sgcca.centroid$a[["metadata"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!duplicated(b[, 1:10]), ]

pvalue <- numeric(ncol(b))
for (col in seq_len(ncol(b))) {
  pvalue[col] <- m(d[col], b[, col])
}
dnames(pvalue) <- names(d)
fdr <- p.adjust(pvalue, "BH")
hist(fdr, main = "Histogram of FDR for metadata", xlab = "FDR")
abline(v = 0.05, col = "red")
fdr[fdr < 0.05]
