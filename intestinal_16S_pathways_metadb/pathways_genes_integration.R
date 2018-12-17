sgcca.centroid <- readRDS("sgcca.RDS")

pathways <- sgcca.centroid$a$RNAseq[, 1]
pathways <- pathways[pathways != 0]

library("reactome.db")
genes_e <- select(reactome.db, keys = names(pathways), keytype = "PATHID", column = "ENTREZID")
pathways_names <- select(reactome.db, keys = names(pathways), keytype = "PATHID", column = "PATHNAME")
pathways_names$PATHNAME <- gsub("(Homo sapiens: )", "", pathways_names$PATHNAME)
library("org.Hs.eg.db")
genes_s <- select(org.Hs.eg.db, keys = genes_e$ENTREZID, keytype = "ENTREZID", column = "SYMBOL")
paths <- merge(genes_e, pathways_names)
paths <- merge(paths, genes_s)
genes_s <- select(org.Hs.eg.db, keys = genes_e$ENTREZID, keytype = "ENTREZID", column = "ENSEMBL")
sgcca.centroid <- readRDS("../intestinal_16S_RNAseq_metadb/sgcca.RDS")
genes <- sgcca.centroid$a$RNAseq[, 1]
names(genes) <- gsub("(ENSG[0-9]+).+", "\\1", names(genes))
genes_se <- unique(genes_s$ENSEMBL)
genes_se <- genes_se[genes_se %in% names(genes)]
genes[genes_se]
(out <- mapIds(org.Hs.eg.db, keys = names(genes[genes_se][genes[genes_se] != 0]), 
       keytype = "ENSEMBL", column = "SYMBOL"))
out <- as.data.frame(out)
write.csv(out, file = "relevant_genes_in_pathway_integration.csv", row.names = TRUE)

micro <- sgcca.centroid$a$`16S`[, 1]
otus <- names(micro)[micro != 0]

cd <- setwd("..")
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

file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(cd)

write.csv(unique(otus_tax_i[otus, ]), file = "microorganisms.csv", na = "", row.names = FALSE)
paths_simple <- paths[, -c(1, 2)]
paths_simple <- paths_simple[order(paths_simple$PATHNAME), ]
write.csv(paths_simple, file = "pathways_genes.csv", row.names = FALSE)
