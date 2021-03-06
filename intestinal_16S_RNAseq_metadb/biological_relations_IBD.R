library("fgsea") # To test the distribution of the pathways
library("reactome.db") # To get the pathways
library("org.Hs.eg.db") # To get the names to entrez from the ensembl
library("data.table")
library("integration")
library("ReactomePA")
library("clusterProfiler")

# Test if the loadings of the genes are from a relevant pathways.
wd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

intestinal <- "intestinal_16S"
stool <- "stools_16S"
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

epithelium <- read.csv("epithelium.csv")
epithelium <- epithelium$Epithelium

setwd(wd)

# Load previous data
sgcca.centroid <- readRDS("IBD.RDS")
STAB <- readRDS("bootstrap_IBD.RDS")

# RNAseq ####
b <- STAB[["RNAseq"]]
d <- sgcca.centroid$a[["RNAseq"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!is.na(b[, 1]), ]

warning(sum(is.na(b[, 1])), " iterations failed.")

pvalue <- numeric(ncol(b))
for (col in seq_len(ncol(b))) {
  pvalue[col] <- two.sided(d[col], b[, col])
}
names(pvalue) <- colnames(b)
fdr <- p.adjust(pvalue, "BH")

# Select those genes that are significant
significant <- names(pvalue)[pvalue < 0.05]
# The fdr results in none significant
significant <- names(pvalue)[pvalue < 0.05]
significant <- gsub("(.*)\\..*", "\\1", significant)

loadings <- sgcca.centroid$a[["RNAseq"]]

# Convert the ids to the entrezIds
ensemblID <- rownames(loadings)
ensemblID <- gsub("(.*)\\..*", "\\1", ensemblID)
rownames(loadings) <- gsub("(.*)\\..*", "\\1", rownames(loadings))
entrezID <- mapIds(
  org.Hs.eg.db, keys = ensemblID, keytype = "ENSEMBL",
  column = "ENTREZID"
)
comp1 <- loadings[, 1]
names(comp1) <- entrezID

epitheliumE <- mapIds(
  org.Hs.eg.db, keys = as.character(epithelium),
  keytype = "SYMBOL", column = "ENTREZID"
)

epitheliumE <- unlist(epitheliumE, use.names = TRUE)
epitheliumE <- epitheliumE[!is.na(epitheliumE)]

# Extract the information of the pathways
genes2Pathways <- as.list(reactome.db::reactomeEXTID2PATHID)
pathways <- unlist(genes2Pathways, use.names = FALSE)
genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
paths2genes <- split(genes, pathways) # List of genes and the gene sets

# Subset to only human pathways
paths2genes <- paths2genes[grep("R-HSA-", names(paths2genes))]

## Compute the hypergeometric/enrichment analysis ####
enrich <- enrichPathway(
  gene = entrezID[significant], pvalueCutoff = 0.05,
  readable = TRUE, universe = unique(entrezID)
)
write.csv(as.data.frame(enrich), file = "RNAseq_enrichment_IBD.csv")

# Store the entrezid
entrezSig <- entrezID[significant]
entrezSig <- entrezSig[!is.na(entrezSig)]
paths2genes[["significant"]] <- entrezSig
paths2genes[["Epithelium"]] <- epitheliumE

## Compute the GSEA for the size effect ####
gseaSizeEffect <- fgsea(paths2genes, comp1, nperm = length(comp1))

# Get the name of the pathway
namesPaths <- select(
  reactome.db, keys = gseaSizeEffect$pathway,
  keytype = "PATHID", columns = "PATHNAME"
)
# Remove the homo sapiens part
namesPaths$PATHNAME <- gsub("Homo sapiens: (.*)", "\\1", namesPaths$PATHNAME)
# Add a column
gseaSizeEffect[, namesPaths := namesPaths$PATHNAME]
# Order the dataframe by size effect
data.table::setorder(gseaSizeEffect, -NES, padj, -size)
# Store the output
fwrite(gseaSizeEffect[padj < 0.05, ], file = "gsea_RNAseq_pathways_IBD.csv")

## 16S ####
b <- STAB[["16S"]]
d <- sgcca.centroid$a[["16S"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!is.na(b[, 1]), ]

pvalue <- numeric(ncol(b))
for (col in seq_len(ncol(b))) {
  pvalue[col] <- two.sided(d[col], b[, col])
}
names(pvalue) <- names(d)
fdr <- p.adjust(pvalue, "BH")

otus <- names(pvalue)[pvalue < 0.05]

## Split the taxonomy
Taxon2Class <- as.list(as.data.frame(otus_tax_i))

grouping <- split(Taxon2Class$Genus, Taxon2Class$Genus)
grouping <- sapply(grouping, names)

term2gene <- data.frame(
  "Gene" = otus_tax_i[, "Genus"],
  "Term" = rownames(otus_tax_i)
)
term2name <- data.frame(
  "Name" = otus_tax_i[, "Genus"],
  "Term" = rownames(otus_tax_i)
)
enrich <- as.data.frame(enricher(
  gene = otus, universe = rownames(otus_tax_i),
  minGSSize = 1, TERM2GENE = term2gene,
  TERM2NAME = term2name
))

write.csv(enrich, file = "Otus_genus_enrichment_IBD.csv")

# GSEA
comp1 <- sgcca.centroid$a[["16S"]][, 1]

gseaSizeEffect <- fgsea(grouping, comp1, nperm = 10000)
data.table::setorder(gseaSizeEffect, -NES, padj, -size)
fwrite(gseaSizeEffect[pval < 0.05], file = "gsea_otus_genus_IBD.csv")
