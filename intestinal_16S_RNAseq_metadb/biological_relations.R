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
all_bootstrap <- readRDS("bootstrap.RDS")
all_sgcca <- readRDS("sgcca.RDS")
write.csv(integration::weights(all_sgcca), file = "RNAseq_weight_all.csv", na = "", row.names = FALSE)
write.csv(integration::weights_otus(all_sgcca, otus_tax_i), file = "16S_weight_all.csv", na = "", row.names = FALSE)
biological_relationships(all_sgcca, all_bootstrap, "all", otus_tax_i, epithelium)

# Controls
controls_bootstrap <- readRDS("bootstrap_Controls.RDS")
controls_sgcca <- readRDS("Controls.RDS")
write.csv(integration::weights(controls_sgcca), file = "RNAseq_weight_controls.csv", na = "", row.names = FALSE)
write.csv(weights_otus(controls_sgcca, otus_tax_i), file = "16S_weight_controls.csv", na = "", row.names = FALSE)
biological_relationships(controls_sgcca, controls_bootstrap, "controls", otus_tax_i, epithelium)

# IBD
IBD_bootstrap <- readRDS("bootstrap_IBD.RDS")
IBD_sgcca <- readRDS("IBD.RDS")
write.csv(integration::weights(IBD_sgcca), file = "RNAseq_weight_IBD.csv", na = "", row.names = FALSE)
write.csv(weights_otus(IBD_sgcca, otus_tax_i), file = "16S_weight_IBD.csv", na = "", row.names = FALSE)
biological_relationships(IBD_sgcca, IBD_bootstrap, "IBD", otus_tax_i, epithelium)

# # T0
# readRDS("bootstrap_IBD_T0.RDS")
# readRDS("IBD_T0.RDS")
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T0.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T0.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T0", otus_tax_i, epithelium)
# 
# # T26
# readRDS("bootstrap_IBD_T26.RDS")
# readRDS("IBD_T26.RDS")
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T26.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T26.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T26", otus_tax_i, epithelium)
# 
# # T52
# readRDS("bootstrap_IBD_T52.RDS")
# readRDS("IBD_T52.RDS")
# write.csv(weights(sgcca.centroid), file = "RNAseq_weight_IBD_T52.csv", na = "", row.names = FALSE)
# write.csv(weights_otus(sgcca.centroid, otus_tax_i), file = "16S_weight_T52.csv", na = "", row.names = FALSE)
# biological_relationships(sgcca.centroid, STAB, "T52", otus_tax_i, epithelium)
