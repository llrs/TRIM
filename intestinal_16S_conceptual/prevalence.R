library("boot")
# Summarize by taxa
library("metagenomeSeq")

wd <- setwd("..")
source("helper_functions.R")
source("helper_prevalence.R")

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
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
setwd(wd)

# Clean the metadata
meta_i <- meta_i_norm(meta_i)

# Match the name of meta and rownames
otus_table_i <- otus_table_i[, match(rownames(meta_i), colnames(otus_table_i))]

MR_i <- newMRexperiment(
  otus_table_i,
  phenoData = AnnotatedDataFrame(meta_i),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)


genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")
species_i <- aggTax(MR_i, lvl = "Species", out = "matrix")
family_i <- aggTax(MR_i, lvl = "Family", out = "matrix")
order_i <- aggTax(MR_i, lvl = "Order", out = "matrix")
class_i <- aggTax(MR_i, lvl = "Class", out = "matrix")
phylum_i <- aggTax(MR_i, lvl = "Phylum", out = "matrix")

# Compare at the otus level the time and response effect
otus <- response_time(otus_table_i, meta_i)
write.csv(otus, "prevalence_otus.csv")
genus_i <- response_time(genus_i, meta_i)
write.csv(genus_i, "prevalence_genus_i.csv")
species_i <- response_time(species_i, meta_i)
write.csv(species_i, "prevalence_species.csv")
family_i <- response_time(family_i, meta_i)
write.csv(family_i, "prevalence_family.csv")
order_i <- response_time(order_i, meta_i)
write.csv(order_i, "prevalence_order.csv")
class_i <- response_time(class_i, meta_i)
write.csv(class_i, "prevalence_class.csv")
phylum_i <- response_time(phylum_i, meta_i)
write.csv(phylum_i, "prevalence_phylum.csv")

# Pairwise comparisons
otus <- comb_prevalence(otus_table_i, meta_i, c("Time", "HSCT_responder"))
write.csv(otus, "otus_pairwise.csv")
species <- comb_prevalence(species_i, meta_i, c("Time", "HSCT_responder"))
write.csv(species, "species_pairwise.csv")
genus <- comb_prevalence(genus_i, meta_i, c("Time", "HSCT_responder"))
write.csv(genus, "genus_pairwise.csv")
family <- comb_prevalence(family_i, meta_i, c("Time", "HSCT_responder"))
write.csv(family, "family_pairwise.csv")
order <- comb_prevalence(order_i, meta_i, c("Time", "HSCT_responder"))
write.csv(order, "order_pairwise.csv")
class <- comb_prevalence(class_i, meta_i, c("Time", "HSCT_responder"))
write.csv(class, "class_pairwise.csv")
phylum <- comb_prevalence(phylum_i, meta_i, c("Time", "HSCT_responder"))
write.csv(phylum, "phylum_pairwise.csv")


# summary:  Not significative, IBD patients are as  likely to have one
# microorganism as the controls.
# Test the prevalence between non controls at times T0, T26, T52 ####
removeControls <- meta_i$IBD == "CD"
removeTimes <- meta_i$Time %in% c("T0", "T26", "T52")
keep <- removeControls & removeTimes
res <- prevalence_tab(
  otus_table_i[, keep], meta_i[keep, ],
  c("Time", "HSCT_responder")
)

# Fisher test
Time_Responder <- prevalence(res$presence, res$absence)

Time_Responder[is.na(Time_Responder)] <- 1
Time_Responder <- p.adjust(Time_Responder, "BH")

if (any(Time_Responder < 0.05)) {
  message("Plots!!")
}

# The presence of microorganisms thorough time is not different in the
# overall patients

# Test if the responders and the non responders behave differently along time
# If the presence of microorganisms thorough time is different in responders
# than in non-responders ####
otus <- comp(otus_table_i, meta_i)
write.csv(otus, "otus_time.csv")
species <- comp(species_i, meta_i)
write.csv(species, "species_time.csv")
genus <- comp(genus_i, meta_i)
write.csv(genus, "genus_time.csv")
family <- comp(family_i, meta_i)
write.csv(family, "family_time.csv")
order <- comp(order_i, meta_i)
write.csv(order, "order_time.csv")
class <- comp(class_i, meta_i)
write.csv(class, "class_time.csv")
phylum <- comp(phylum_i, meta_i)
write.csv(phylum, "phylum_time.csv")


# Bootstrapp while controlling by Time ####

# # bootstrapping with 1000 replications
# results <- boot(data = t(sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005),
#                 statistic = ratio,
#                 R = 1000,
#                 meta = meta_i[removeControls, ], columns = "Time",
#                 parallel = "multicore", ncpus = 4)
#
# # view results
# results
# plot(results)
#
# # get 95% confidence interval
# boot.ci(results, type="bca")

dev.off()
