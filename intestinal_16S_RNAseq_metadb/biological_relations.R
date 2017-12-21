# Test if the loadings of the genes are from a relevant pathways. 

# Load the helper file
source("helper_functions.R")

library("fgsea") # To test the distribution of the pathways
library("reactome.db") # To get the pathways
library("org.Hs.eg.db") # To get the names to entrez from the ensembl

# Load previous data
load("sgcca.RData")

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

# Compute the GSEA
f <- fgsea(paths2genes, comp1, nperm = length(comp1))
signPathways <- f$pathway[f$padj < 0.05]
namesPaths <- select(reactome.db, keys = signPathways, 
                     keytype = "PATHID", columns = "PATHNAME")
out <- cbind(f[f$padj < 0.05, ], namesPaths = namesPaths$PATHNAME)
data.table::setorder(out, -NES, padj, -size)
write.csv(out, file = "relevant_pathways.csv")
