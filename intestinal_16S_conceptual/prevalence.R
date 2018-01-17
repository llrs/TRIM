library("boot")
# Summarize by taxa
library("metagenomeSeq")

wd <- setwd("..")
source("helper_functions.R")
source("helper_prevalence.R")

intestinal <- "intestinal_16S"
# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
setwd(wd)


# Clean the metadata
meta_i <- meta_i[, apply(meta_i, 2, function(x){length(unique(x)) != 1})]
meta_i$ID <- meta_i$Patient_ID
meta_i$ID[meta_i$Patient_ID %in% c("15", "23")] <- "15/23"
meta_i$ID[meta_i$Patient_ID %in% c("33", "36")] <- "33/36"
meta_i$ID[meta_i$Patient_ID %in% c("29", "35")] <- "29/35"
meta_i$ID <- as.factor(meta_i$ID)

# There is a mislabeling on those tubes, we don't know which is which
meta_i$CD_Aftected_area[meta_i$Sample_Code == "22_T52_T_DM_III"] <- NA

meta_i$HSCT_responder[(meta_i$ID %in% c("38", "40", "41"))] <- NA

# Match the name of meta and rownames
otus_table_i <- otus_table_i[, match(rownames(meta_i), colnames(otus_table_i))]

MR_i <- newMRexperiment(otus_table_i, 
                        phenoData = AnnotatedDataFrame(meta_i),
                        featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i)))
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")


# Test the prevalence between controls and non controls ####
res <- prevalence_tab(otus_table_i, meta_i, "IBD")

# Fisher test
IBD <- prevalence(res$presence, res$absence)
  
IBD[is.na(IBD)] <- 1
IBD <- p.adjust(IBD, "BH")

summary(IBD < 0.05)

pdf(paste0("Figures/", today, "_ratios_plots.pdf"))

# Differences between responders and non responders at time 0 ####
removeControls <- meta_i$IBD  != "CONTROL"
removeTimes <- meta_i$Time  == "T0"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_i[, keep], meta_i[keep, ], "HSCT_responder")

# Fisher test
T0 <- prevalence(res$presence, res$absence)
T0[is.na(T0)] <- 1
T0 <- p.adjust(T0, "BH")

summary(T0 < 0.05)

# Differences between responders and non responders at time T26 ####
removeControls <- meta_i$IBD  != "CONTROL"
removeTimes <- meta_i$Time  == "T26"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_i[, keep], meta_i[keep, ], "HSCT_responder")

# Fisher test
T26 <- prevalence(res$presence, res$absence)
T26[is.na(T26)] <- 1
T26 <- p.adjust(T26, "BH")

summary(T26 < 0.05)

# Differences between responders and non responders at time T52 ####
removeControls <- meta_i$IBD  != "CONTROL"
removeTimes <- meta_i$Time  == "T52"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_i[, keep], meta_i[keep, ], "HSCT_responder")

# Fisher test
T52 <- prevalence(res$presence, res$absence)
T52[is.na(T52)] <- 1
T52 <- p.adjust(T52, "BH")

summary(T52 < 0.05)

# summary:  Not significative, IBD patients are as  likely to have one 
# microorganism as the controls.

# Test the prevalence between non controls at times T0, T26, T52 ####
removeControls <- meta_i$IBD  == "CD"
removeTimes <- meta_i$Time %in% c("T0", "T26", "T52")
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_i[, keep], meta_i[keep, ], 
                      c("Time", "HSCT_responder"))

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

species_i <- aggTax(MR_i, lvl = "Species", out = "matrix")
species <- comp(species_i, meta_i)
unique(otus_tax_i[otus_tax_i[, "Species"] %in% names(which(species$Responders < 0.05)), 
                  c("Genus", "Species")])

genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")
genus <- comp(genus_i, meta_i)
which(genus$Responders < 0.05)

family_i <- aggTax(MR_i, lvl = "Family", out = "matrix")
family <- comp(family_i, meta_i)
which(family$Responders < 0.05)

order_i <- aggTax(MR_i, lvl = "Order", out = "matrix")
order <- comp(order_i, meta_i)

class_i <- aggTax(MR_i, lvl = "Class", out = "matrix")
class <- comp(class_i, meta_i)

phylum_i <- aggTax(MR_i, lvl = "Phylum", out = "matrix")
phylum <- comp(phylum_i, meta_i)

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