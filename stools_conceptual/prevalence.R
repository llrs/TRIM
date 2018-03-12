library("boot")
# Summarize by taxa
library("metagenomeSeq")

wd <- setwd("..")
library("integration")

stool <- "stools_16S"
# Read the stools OTUs
otus_table_s <- read.delim(
  file.path(stool, "OTUs-Table-refined-stools.tab"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_s <- otus_table_s[, ncol(otus_table_s)]
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

# Extract the taxonomy and format it properly
otus_tax_s <- taxonomy(tax_s, rownames(otus_table_s))

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(
  file_meta_s, check.names = FALSE, row.names = 1,
  stringsAsFactors = FALSE
)
setwd(wd)

# Reorder columnames to match rownames of metadata
otus_table_s <- otus_table_s[, match(rownames(meta_s), colnames(otus_table_s))]

MR_s <- newMRexperiment(
  otus_table_s,
  phenoData = AnnotatedDataFrame(meta_s),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_s))
)

genus_i <- aggTax(MR_s, lvl = "Genus", out = "matrix")
species_i <- aggTax(MR_s, lvl = "Species", out = "matrix")
family_i <- aggTax(MR_s, lvl = "Family", out = "matrix")
order_i <- aggTax(MR_s, lvl = "Order", out = "matrix")
class_i <- aggTax(MR_s, lvl = "Class", out = "matrix")
phylum_i <- aggTax(MR_s, lvl = "Phylum", out = "matrix")

pdf(paste0("Figures/", today, "_ratios_plots.pdf"))

# Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x) {
  length(unique(x)) != 1
})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)



## Time ####
otus <- comb_prevalence(otus_table_s, meta_s, c("Time"))
write.csv(otus, "prevalence_time_otus.csv")
genus <- comb_prevalence(genus_i, meta_s, c("Time"))
write.csv(genus, "prevalence_time_genus.csv")
species <- comb_prevalence(species_i, meta_s, c("Time"))
write.csv(species, "prevalence_time_species.csv")
family <- comb_prevalence(family_i, meta_s, c("Time"))
write.csv(family, "prevalence_time_family.csv")
order <- comb_prevalence(order_i, meta_s, c("Time"))
write.csv(order, "prevalence_time_order.csv")
class <- comb_prevalence(class_i, meta_s, c("Time"))
write.csv(class, "prevalence_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta_s, c("Time"))
write.csv(phylum, "prevalence_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_table_s, meta_s, c("HSCT_responder"))
write.csv(otus, "prevalence_otus.csv")
genus <- comb_prevalence(genus_i, meta_s, c("HSCT_responder"))
write.csv(genus, "prevalence_genus_i.csv")
species <- comb_prevalence(species_i, meta_s, c("HSCT_responder"))
write.csv(species, "prevalence_species.csv")
family <- comb_prevalence(family_i, meta_s, c("HSCT_responder"))
write.csv(family, "prevalence_family.csv")
order <- comb_prevalence(order_i, meta_s, c("HSCT_responder"))
write.csv(order, "prevalence_order.csv")
class <- comb_prevalence(class_i, meta_s, c("HSCT_responder"))
write.csv(class, "prevalence_class.csv")
phylum <- comb_prevalence(phylum_i, meta_s, c("HSCT_responder"))
write.csv(phylum, "prevalence_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_table_s, meta_s, c("HSCT_responder", "Time"))
write.csv(otus, "prevalence_response_time_otus.csv")
genus <- comb_prevalence(genus_i, meta_s, c("HSCT_responder", "Time"))
write.csv(genus, "prevalence_response_time_genus_i.csv")
species <- comb_prevalence(species_i, meta_s, c("HSCT_responder", "Time"))
write.csv(species, "prevalence_response_time_species.csv")
family <- comb_prevalence(family_i, meta_s, c("HSCT_responder", "Time"))
write.csv(family, "prevalence_response_time_family.csv")
order <- comb_prevalence(order_i, meta_s, c("HSCT_responder", "Time"))
write.csv(order, "prevalence_response_time_order.csv")
class <- comb_prevalence(class_i, meta_s, c("HSCT_responder", "Time"))
write.csv(class, "prevalence_response_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta_s, c("HSCT_responder", "Time"))
write.csv(phylum, "prevalence_response_time_phylum.csv")

# No need to test for different locations :D

dev.off()
