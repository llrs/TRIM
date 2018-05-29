# Test if the loadings of the genes are from a relevant pathways.
today <- format(Sys.time(), "%Y%m%d")
library("integration")
wd <- setwd("..")

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

# Do the biological analysis of all the samples
load("bootstrap.RData", verbose = TRUE)
load("sgcca.RData", verbose = TRUE)
write.csv(integration::weights(sgcca.centroid), file = "RNAseq_weight_all.csv", na = "", row.names = FALSE)
write.csv(integration::weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_all.csv", na = "", row.names = FALSE)
biological_relationships(sgcca.centroid, STAB, "all", otus_tax_i, epithelium)

# Controls
load("bootstrap_Controls.RData", verbose = TRUE)
load("Controls.RData", verbose = TRUE)
write.csv(integration::weights(sgcca.centroid), file = "RNAseq_weight_controls.csv", na = "", row.names = FALSE)
write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_controls.csv", na = "", row.names = FALSE)
biological_relationships(sgcca.centroid, STAB, "controls", otus_tax_i, epithelium)

# IBD
load("bootstrap_IBD.RData", verbose = TRUE)
load("IBD.RData", verbose = TRUE)
write.csv(integration::weights(sgcca.centroid), file = "RNAseq_weight_IBD.csv", na = "", row.names = FALSE)
write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_IBD.csv", na = "", row.names = FALSE)
biological_relationships(sgcca.centroid, STAB, "IBD", otus_tax_i, epithelium)

# # T0
# load("bootstrap_IBD_T0.RData", verbose = TRUE)
# load("IBD_T0.RData", verbose = TRUE)
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T0.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T0.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T0", otus_tax_i, epithelium)
# 
# # T26
# load("bootstrap_IBD_T26.RData", verbose = TRUE)
# load("IBD_T26.RData", verbose = TRUE)
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T26.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T26.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T26", otus_tax_i, epithelium)
# 
# # T52
# load("bootstrap_IBD_T52.RData", verbose = TRUE)
# load("IBD_T52.RData", verbose = TRUE)
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T52.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T52.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T52", otus_tax_i, epithelium)
