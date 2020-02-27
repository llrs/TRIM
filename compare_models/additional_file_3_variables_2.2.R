library("org.Hs.eg.db")
library("integration")

intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"

folder1 <- "intestinal_16S_RNAseq_metadb"
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)

tax_i <- otus_table_i[, ncol(otus_table_i)]
ta <- integration::taxonomy(tax_i, rownames(otus_table_i))

# * Export weight of model 2.2 ####
r <- apply(model2.2$a$RNAseq, 1, function(x){any(x != 0)})
ids <- mapIds(org.Hs.eg.db, keys = trimVer(names(r)[r]), column = "SYMBOL", 
              keytype = "ENSEMBL")
out <- cbind(Symbol = ids, model2.2$a$RNAseq[r, ])
write.csv(out, file = "data_out/m2.2_RNA_weights.csv")
r <- apply(model2.2$a$`16S`, 1, function(x){any(x != 0)})
out <- cbind(ta[names(r)[r], ], model2.2$a$`16S`[r, ])
write.csv(out, file = "data_out/m2.2_16S_weights.csv")
