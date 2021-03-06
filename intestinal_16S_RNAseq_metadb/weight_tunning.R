cd <- setwd("..")
library("fgsea")
library("org.Hs.eg.db")
library("reactome.db")
library("ggforce")
library("metagenomeSeq")
library("integration")
library("ggplot2")
library("clusterProfiler")
library("KEGGREST")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

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

# Correct the swapped samples and match metadata
expr <- norm_expr_colnames(expr)
# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

expr_norm <- norm_RNAseq(expr)
expr <- filter_RNAseq(expr_norm)
otus_table_i <- norm_otus(otus_table_i)

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them
metadb <- meta_r
keepCol <- sapply(metadb, is.factor)
nam <- c(
  "Exact_location", # Segment of the sample
  # superseeded by SESCD 
  # "Active_area", # Health stage of the sample
  # "IBD", # Disease or control
  "AGE_SAMPLE", # Age
  # "diagTime", # Time with disease
  # Not really needed induced by diagTime and age sample
  "AgeDiag", # Age at which the disease was diagnositcated
  "Transplant", # Stage of the treatment
  "ID", # Patient
  # "SESCD_local", # Clinical score
  "Treatment", # Further complications
  "Surgery", # Up to surgery?
  "SEX" # Male/female
) 

o <- sapply(nam[-c(2, 3)], allComb, data = droplevels(meta_r))
ou <- sapply(meta_r[, nam[-c(2, 3)]], function(x){length(unique(x))})
o[ou == 2] <- lapply(o[ou == 2], function(x){x[, 1]})
o2 <- do.call(base::cbind, o)
o2[o2 == TRUE] <- 1
o2[o2 == FALSE] <- 0
out <- cbind(o2, meta_r[, nam[c(2, 3)]])
out[55, ][is.na(out[55, ])] <- c(1, 1, 0, 1, 1, 1)
out[54, ][is.na(out[54, ])] <- c(0, 0, 1, 0, 0, 0)
out[, 48][is.na(out[, 48])] <- out[, 47][is.na(out[, 48])]

metadb <- model_RGCCA(meta_r, nam) 
# Set metadb with a sigle variable with several options
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), "metadata" = metadb)
A <- clean_unvariable(A)
saveRDS(A, file = "TRIM.RDS")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "metadata", 1)
C <- subSymm(C, "RNAseq", "metadata", 1)

weights <- seq(0.05, 1, by = 0.05)
# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.25670333, 0, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

results <- vector(length = 20*20, mode = "list")
names(results) <- apply(expand.grid(weights, weights), 1, paste, collapse = "_")

for (i in weights) {
  C <- subSymm(C, "RNAseq", "metadata", i)
  for (j in weights) {
    C <- subSymm(C, "16S", "metadata", j)
    # print(C)
    sgcca.centroid <- sgcca(
      A, C, c1 = shrinkage,
      ncomp = ncomp,
      scheme = "centroid",
      scale = TRUE,
      verbose = FALSE
    )
    sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))
    results[[paste(i, j, sep = "_")]] <- sgcca.centroid
    
  }
}

saveRDS(results, file = "weights_design.RDS")
nam_genes <- rownames(results[[1]]$a$RNAseq)
nam_genes <- gsub("\\..+", "", nam_genes)
nam_entrez <- mapIds(org.Hs.eg.db, nam_genes, keytype = "ENSEMBL", "ENTREZID")
genes2Pathways <- as.list(reactome.db::reactomeEXTID2PATHID)
pathways <- unlist(genes2Pathways, use.names = FALSE)
genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
paths2genes <- split(genes, pathways) # List of genes and the gene sets

# Subset to only human pathways
paths2genes <- paths2genes[grep("R-HSA-", names(paths2genes))]

o <- vapply(results, function(x){
  genes <- x$a$RNAseq[, 1]
  names(genes) <- nam_entrez
  nrow(enrichKEGG(gene         = names(genes)[genes != 0],
                  organism     = 'hsa',
                  pvalueCutoff = 0.05))
}, numeric(1L))


output <- t(simplify2array(strsplit(names(o), "_")))
output <- cbind(output, enrichment = o)
colnames(output) <- c("RNAseq", "Micro", "Enrichment")
output <- apply(output, 2, as.numeric)
output <- as.data.frame(output)
output$genes <- vapply(results, function(x){
  genes <- x$a$RNAseq[, 1]
  sum(genes != 0, na.rm = TRUE)
  }, numeric(1L))

ggplot(output, aes(x = RNAseq, y = Micro)) + 
  geom_count(aes(size = Enrichment, color = Enrichment)) +
  theme_bw() +
  ggtitle("Values in design matrix")
ggplot(output, aes(x = RNAseq, y = Micro)) + 
  geom_contour(aes(z = Enrichment, col = Enrichment)) +
  theme_bw() +
  ggtitle("Values in design matrix")

hsa_path_eg  <- keggLink("pathway", "hsa")
hsa_eg <- gsub("hsa:", "", names(hsa_path_eg))
hsa_en <- mapIds(org.Hs.eg.db, hsa_eg, keytype = "ENTREZID", "ENSEMBL")
paths2genes <- split(hsa_en, hsa_path_eg)
paths2genes2 <- lapply(paths2genes, unique)


o <- vapply(results, function(x){
  genes <- x$a$RNAseq[, 1]
  # names(genes) <- nam_entrez
  out <- fgsea(paths2genes2, genes[genes != 0], nperm = 5000)
  sum(out$padj < 0.05)
}, numeric(1L))


output <- t(simplify2array(strsplit(names(o), "_")))
output <- cbind(output, enrichment = o)
colnames(output) <- c("RNAseq", "Micro", "Enrichment")
output <- apply(output, 2, as.numeric)
output <- as.data.frame(output)
ggplot(output, aes(x = RNAseq, y = Micro)) + 
  geom_count(aes(size = Enrichment, color = Enrichment)) +
  theme_bw() +
  ggtitle("Values in design matrix")
ggplot(output, aes(x = RNAseq, y = Micro)) + 
  geom_contour(aes(z = Enrichment, col = Enrichment)) +
  theme_bw() +
  ggtitle("Values in design matrix")
