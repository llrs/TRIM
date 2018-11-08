library("integration")
library("fgsea")
library("ggforce")
library("RGCCA")
library("reactome.db")
library("org.Hs.eg.db")

today <- format(Sys.time(), "%Y%m%d")

# Save
otus_table_i <- readRDS("otus_table.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS( "meta.RDS")

epithelium <- read.csv("../epithelium.csv")
epithelium <- epithelium$Epithelium

# Translate into Entrez
epithelium <- epitheliumE(epithelium)

# Do the biological analysis of all the samples
all_bootstrap <- readRDS("bootstrap.RDS")
all_sgcca <- readRDS("sgcca.RDS")
write.csv(integration::weights(all_sgcca), file = "RNAseq_weight_all.csv", na = "", row.names = FALSE)
write.csv(integration::weights_otus(all_sgcca, otus_tax_i), file = "16S_weight_all.csv", na = "", row.names = FALSE)
biological_relationships(all_sgcca, all_bootstrap$STAB, "all", otus_tax_i, epithelium, today)

# Controls
controls_bootstrap <- readRDS("bootstrap_Controls.RDS")
controls_sgcca <- readRDS("Controls.RDS")
write.csv(integration::weights(controls_sgcca), file = "RNAseq_weight_controls.csv", na = "", row.names = FALSE)
write.csv(weights_otus(controls_sgcca, otus_tax_i), file = "16S_weight_controls.csv", na = "", row.names = FALSE)
biological_relationships(controls_sgcca, controls_bootstrap$STAB, "controls", otus_tax_i, epithelium, today)

# IBD
IBD_bootstrap <- readRDS("bootstrap_IBD.RDS")
IBD_sgcca <- readRDS("IBD.RDS")
write.csv(integration::weights(IBD_sgcca), file = "RNAseq_weight_IBD.csv", na = "", row.names = FALSE)
write.csv(weights_otus(IBD_sgcca, otus_tax_i), file = "16S_weight_IBD.csv", na = "", row.names = FALSE)
biological_relationships(IBD_sgcca, IBD_bootstrap$STAB, "IBD", otus_tax_i, epithelium, today)
