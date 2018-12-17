library("boot")
# Summarize by taxa
library("metagenomeSeq")

wd <- setwd("..")
today <- format(Sys.time(), "%Y%m%d")
library("integration")

intestinal <- "intestinal_16S"
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

# Read the metadata for each type of sample
rna <- "intestinal_RNAseq"
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(wd)

# Clean the metadata
meta_r <- meta_r_norm(meta_r)

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Subset and reorder
meta_r <- meta_r[match(colnames(otus_table_i), meta_r$Seq_code_uDNA), ]

# Subset only the involved
involve <- meta_r$Involved_Healthy %in% "INVOLVED" | meta_r$IBD %in% "CONTROL"

otus_table_i <- otus_table_i[, involve]
otus_tax_i <- otus_tax_i[rowSums(otus_table_i) != 0, ]
otus_table_i <- otus_table_i[rowSums(otus_table_i) != 0, ]
meta_r <- meta_r[involve, ]

rownames(meta_r) <- meta_r$Seq_code_uDNA

MR_i <- newMRexperiment(
  otus_table_i,
  phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)


genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")
species_i <- aggTax(MR_i, lvl = "Species", out = "matrix")
family_i <- aggTax(MR_i, lvl = "Family", out = "matrix")
order_i <- aggTax(MR_i, lvl = "Order", out = "matrix")
class_i <- aggTax(MR_i, lvl = "Class", out = "matrix")
phylum_i <- aggTax(MR_i, lvl = "Phylum", out = "matrix")

# Old ####
# # Compare at the otus level the time and response effect 
# # otus <- response_time(otus_table_i, meta_r)
# # write.csv(otus, "prevalence_otus.csv")
# # genus <- response_time(genus_i, meta_r)
# # write.csv(genus, "prevalence_genus_i.csv")
# # species <- response_time(species_i, meta_r)
# # write.csv(species, "prevalence_species.csv")
# # family <- response_time(family_i, meta_r)
# # write.csv(family, "prevalence_family.csv")
# # order <- response_time(order_i, meta_r)
# # write.csv(order, "prevalence_order.csv")
# # class <- response_time(class_i, meta_r)
# # write.csv(class, "prevalence_class.csv")
# phylum <- response_time(phylum_i, meta_r)
# write.csv(phylum, "prevalence_phylum.csv")
# 
# # Pairwise comparisons
# otus <- comb_prevalence(otus_table_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(otus, "otus_pairwise.csv")
# species <- comb_prevalence(species_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(species, "species_pairwise.csv")
# genus <- comb_prevalence(genus_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(genus, "genus_pairwise.csv")
# family <- comb_prevalence(family_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(family, "family_pairwise.csv")
# order <- comb_prevalence(order_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(order, "order_pairwise.csv")
# class <- comb_prevalence(class_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(class, "class_pairwise.csv")
# phylum <- comb_prevalence(phylum_i, meta_r, c("Time", "HSCT_responder"))
# write.csv(phylum, "phylum_pairwise.csv")
# 
# 
# # summary:  Not significative, IBD patients are as  likely to have one
# # microorganism as the controls.
# # Test the prevalence between non controls at times T0, T26, T52 #
# removeControls <- meta_r$IBD == "CD"
# removeTimes <- meta_r$Time %in% c("T0", "T26", "T52")
# keep <- removeControls & removeTimes
# Time_Responder <- comb_prevalence(otus_table_i[, keep], meta_r[keep, ],
#   c("Time", "HSCT_responder"))
# 
# if (any(Time_Responder < 0.05)) {
#   message("Plots!!")
# }
# 
# # The presence of microorganisms thorough time is not different in the
# # overall patients
# 
# # Test if the responders and the non responders behave differently along time
# # If the presence of microorganisms thorough time is different in responders
# # than in non-responders #
# otus <- comp(otus_table_i, meta_r)
# write.csv(otus, "otus_time.csv")
# species <- comp(species_i, meta_r)
# write.csv(species, "species_time.csv")
# genus <- comp(genus_i, meta_r)
# write.csv(genus, "genus_time.csv")
# family <- comp(family_i, meta_r)
# write.csv(family, "family_time.csv")
# order <- comp(order_i, meta_r)
# write.csv(order, "order_time.csv")
# class <- comp(class_i, meta_r)
# write.csv(class, "class_time.csv")
# phylum <- comp(phylum_i, meta_r)
# write.csv(phylum, "phylum_time.csv")



## Time ####
otus <- comb_prevalence(otus_table_i, meta_r, c("Time"))
write.csv(otus, "prevalence_time_otus.csv")
genus <- comb_prevalence(genus_i, meta_r, c("Time"))
write.csv(genus, "prevalence_time_genus.csv")
species <- comb_prevalence(species_i, meta_r, c("Time"))
write.csv(species, "prevalence_time_species.csv")
family <- comb_prevalence(family_i, meta_r, c("Time"))
write.csv(family, "prevalence_time_family.csv")
order <- comb_prevalence(order_i, meta_r, c("Time"))
write.csv(order, "prevalence_time_order.csv")
class <- comb_prevalence(class_i, meta_r, c("Time"))
write.csv(class, "prevalence_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta_r, c("Time"))
write.csv(phylum, "prevalence_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_table_i, meta_r, c("HSCT_responder"))
write.csv(otus, "prevalence_otus.csv")
genus <- comb_prevalence(genus_i, meta_r, c("HSCT_responder"))
write.csv(genus, "prevalence_genus_i.csv")
species <- comb_prevalence(species_i, meta_r, c("HSCT_responder"))
write.csv(species, "prevalence_species.csv")
family <- comb_prevalence(family_i, meta_r, c("HSCT_responder"))
write.csv(family, "prevalence_family.csv")
order <- comb_prevalence(order_i, meta_r, c("HSCT_responder"))
write.csv(order, "prevalence_order.csv")
class <- comb_prevalence(class_i, meta_r, c("HSCT_responder"))
write.csv(class, "prevalence_class.csv")
phylum <- comb_prevalence(phylum_i, meta_r, c("HSCT_responder"))
write.csv(phylum, "prevalence_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_table_i, meta_r, c("HSCT_responder", "Time"))
write.csv(otus, "prevalence_response_time_otus.csv")
genus <- comb_prevalence(genus_i, meta_r, c("HSCT_responder", "Time"))
write.csv(genus, "prevalence_response_time_genus_i.csv")
species <- comb_prevalence(species_i, meta_r, c("HSCT_responder", "Time"))
write.csv(species, "prevalence_response_time_species.csv")
family <- comb_prevalence(family_i, meta_r, c("HSCT_responder", "Time"))
write.csv(family, "prevalence_response_time_family.csv")
order <- comb_prevalence(order_i, meta_r, c("HSCT_responder", "Time"))
write.csv(order, "prevalence_response_time_order.csv")
class <- comb_prevalence(class_i, meta_r, c("HSCT_responder", "Time"))
write.csv(class, "prevalence_response_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta_r, c("HSCT_responder", "Time"))
write.csv(phylum, "prevalence_response_time_phylum.csv")




# COLON ######

colon <- !meta_r$CD_Aftected_area %in% "ILEUM"

## Time ####
otus <- comb_prevalence(otus_table_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(otus, "prevalence_colon_time_otus.csv")
genus <- comb_prevalence(genus_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(genus, "prevalence_colon_time_genus.csv")
species <- comb_prevalence(species_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(species, "prevalence_colon_time_species.csv")
family <- comb_prevalence(family_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(family, "prevalence_colon_time_family.csv")
order <- comb_prevalence(order_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(order, "prevalence_colon_time_order.csv")
class <- comb_prevalence(class_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(class, "prevalence_colon_time_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], droplevels(meta_r[colon, ]), c("Time"))
write.csv(phylum, "prevalence_colon_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_table_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(otus, "prevalence_colon_otus.csv")
genus <- comb_prevalence(genus_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(genus, "prevalence_colon_genus_i.csv")
species <- comb_prevalence(species_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(species, "prevalence_colon_species.csv")
family <- comb_prevalence(family_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(family, "prevalence_colon_family.csv")
order <- comb_prevalence(order_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(order, "prevalence_colon_order.csv")
class <- comb_prevalence(class_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(class, "prevalence_colon_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder"))
write.csv(phylum, "prevalence_colon_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_table_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(otus, "prevalence_colon_response_time_otus.csv")
genus <- comb_prevalence(genus_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(genus, "prevalence_colon_response_time_genus_i.csv")
species <- comb_prevalence(species_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(species, "prevalence_colon_response_time_species.csv")
family <- comb_prevalence(family_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(family, "prevalence_colon_response_time_family.csv")
order <- comb_prevalence(order_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(order, "prevalence_colon_response_time_order.csv")
class <- comb_prevalence(class_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(class, "prevalence_colon_response_time_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], droplevels(meta_r[colon, ]), c("HSCT_responder", "Time"))
write.csv(phylum, "prevalence_colon_response_time_phylum.csv")


# ILEUM ####

ileum <- meta_r$CD_Aftected_area %in% "ILEUM"
  
## Time ####
otus <- comb_prevalence(otus_table_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(otus, "prevalence_ileum_time_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(genus, "prevalence_ileum_time_genus.csv")
species <- comb_prevalence(species_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(species, "prevalence_ileum_time_species.csv")
family <- comb_prevalence(family_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(family, "prevalence_ileum_time_family.csv")
order <- comb_prevalence(order_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(order, "prevalence_ileum_time_order.csv")
class <- comb_prevalence(class_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(class, "prevalence_ileum_time_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], droplevels(meta_r[ileum, ]), c("Time"))
write.csv(phylum, "prevalence_ileum_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_table_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(otus, "prevalence_ileum_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(genus, "prevalence_ileum_genus_i.csv")
species <- comb_prevalence(species_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(species, "prevalence_ileum_species.csv")
family <- comb_prevalence(family_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(family, "prevalence_ileum_family.csv")
order <- comb_prevalence(order_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(order, "prevalence_ileum_order.csv")
class <- comb_prevalence(class_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(class, "prevalence_ileum_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder"))
write.csv(phylum, "prevalence_ileum_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_table_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(otus, "prevalence_ileum_response_time_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(genus, "prevalence_ileum_response_time_genus_i.csv")
species <- comb_prevalence(species_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(species, "prevalence_ileum_response_time_species.csv")
family <- comb_prevalence(family_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(family, "prevalence_ileum_response_time_family.csv")
order <- comb_prevalence(order_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(order, "prevalence_ileum_response_time_order.csv")
class <- comb_prevalence(class_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(class, "prevalence_ileum_response_time_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], droplevels(meta_r[ileum, ]), c("HSCT_responder", "Time"))
write.csv(phylum, "prevalence_ileum_response_time_phylum.csv")




# Bootstrapp while controlling by Time ####

# # bootstrapping with 1000 replications
# results <- boot(data = t(sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005),
#                 statistic = ratio,
#                 R = 1000,
#                 meta = meta_r[removeControls, ], columns = "Time",
#                 parallel = "multicore", ncpus = 4)
#
# # view results
# results
# plot(results)
#
# # get 95% confidence interval
# boot.ci(results, type="bca")

# dev.off()
