library("ggforce")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")


meta_r$Location <- ifelse(meta_r$Exact_location != "ILEUM", "COLON", "ILEUM")
levels(as.factor(meta_r$Location))
levels(as.factor(meta_r$IBD))
keep <- allComb(meta_r, c("Location", "IBD"))

