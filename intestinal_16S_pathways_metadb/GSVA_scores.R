cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
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

# Load the input data
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

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

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes = TRUE, log = TRUE)

# Filter expression
expr <- norm_RNAseq(expr_norm)

otus_table_i <- otus_table_i[rowSums(otus_table_i) != 0, ]

# Normalize OTUS
library("metagenomeSeq")
MR_i <- newMRexperiment(
  otus_table_i, 
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i[rownames(otus_table_i), ]))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus_table_i <- MRcounts(MR_i, norm = TRUE, log = TRUE)

library("GSVA")
library("reactome.db")
library("GSEABase")
library("org.Hs.eg.db")
paths2genes <- as.list(reactomePATHID2EXTID)
paths2genes <- paths2genes[grep("R-HSA-", names(paths2genes))]
gsl <- sapply(names(paths2genes), function(x){
  try({
    ids <- mapIds(org.Hs.eg.db, keys = unique(paths2genes[[x]]), 
                  keytype = "ENTREZID", column = "ENSEMBL")
    GeneSet(unique(ids), setName = x)}
    , silent = TRUE)
})
l_class <- sapply(gsl, function(x){is(x, "GeneSet")})
gsl <- gsl[l_class]
gsc <- GeneSetCollection(gsl)
rownames(expr) <- gsub("(.*)\\..*", "\\1", rownames(expr))
gsv <- gsva(expr, gset.idx.list = gsc)
saveRDS(gsv, file = "GSV.RDS")

# 
# su <- apply(meta_r[, c("ID", "Exact_location", "AGE_SAMPLE", "diagTime", "SESCD_local", "SEX")], 2, is.na)
# pos <- which(rowSums(su) == 1, arr.ind = TRUE)
# 
# meta_r$SEX <- as.factor(meta_r$SEX)
# 
# library("vegan") # For all the matrice
# all_RNA <- adonis(as.data.frame(t(gsv))[-pos, ] ~ Exact_location * IBD*ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[-pos, ], method = "euclidean")
# 
# # CONTROLS
# contr <- setdiff(which(meta_r$IBD == "CONTROL"), pos)
# 
# c_RNA <- adonis(as.data.frame(t(gsv))[contr, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local + SEX, meta_r[contr, ], method = "euclidean")
# 
# 
# # IBD
# ibd <- setdiff(which(meta_r$IBD != "CONTROL"), pos)
# 
# i_RNA <- adonis(as.data.frame(t(gsv))[ibd, ] ~ Exact_location * ID + AGE_SAMPLE + diagTime + SESCD_local, meta_r[ibd, ], method = "euclidean")
