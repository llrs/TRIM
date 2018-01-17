library("boot")
# Summarize by taxa
library("metagenomeSeq")

wd <- setwd("..")
source("helper_functions.R")
source("helper_prevalence.R")

stool <- "stools_16S"
# Read the stools OTUs
otus_table_s <- read.delim(file.path(stool, "OTUs-Table-refined-stools.tab"), 
                           stringsAsFactors = FALSE, row.names = 1,
                           check.names = FALSE)
tax_s <- otus_table_s[, ncol(otus_table_s)]
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

# Extract the taxonomy and format it properly
otus_tax_s <- taxonomy(tax_s, rownames(otus_table_s))

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(file_meta_s, check.names = FALSE, row.names = 1, 
                     stringsAsFactors = FALSE)
setwd(wd)

# Reorder columnames to match rownames of metadata
otus_table_s <- otus_table_s[, match(rownames(meta_s), colnames(otus_table_s))]

MR_s <- newMRexperiment(otus_table_s, 
                        phenoData = AnnotatedDataFrame(meta_s),
                        featureData = AnnotatedDataFrame(as.data.frame(otus_tax_s)))
genus_s <- aggTax(MR_s, lvl = "Species", out = "matrix")

pdf(paste0("Figures/", today, "_ratios_plots.pdf"))

# Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x){length(unique(x)) != 1})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)

# Test the prevalence between controls and non controls ####
res <- prevalence_tab(otus_table_s, meta_s, "Controls")

# Fisher test
Disease <- prevalence(res$presence, res$absence)
Disease[is.na(Disease)] <- 1
Disease <- p.adjust(Disease, "BH")

# summary:  Not significative, IBD patients are as  likely to have one 
# microorganism as the controls.

# Differences between responders and non responders at time 0
removeControls <- meta_s$Controls  == "NO"
removeTimes <- meta_s$Time  == "T0"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_s[, keep], meta_s[keep, ], "HSCT_responder")

# Fisher test
T0 <- prevalence(res$presence, res$absence)
T0[is.na(T0)] <- 1
T0 <- p.adjust(T0, "BH")

any(T0 < 0.05)

# Differences between responders and non responders at time T26 ##
removeControls <- meta_s$Controls  == "NO"
removeTimes <- meta_s$Time  == "T26"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_s[, keep], meta_s[keep, ], "HSCT_responder")

# Fisher test
T26 <- prevalence(res$presence, res$absence)
T26[is.na(T26)] <- 1
T26 <- p.adjust(T26, "BH")

any(T26 < 0.05)

# Differences between responders and non responders at time T52 ##
removeControls <- meta_s$Controls  == "NO"
removeTimes <- meta_s$Time  == "T52"
keep <- removeControls & removeTimes
res <- prevalence_tab(otus_table_s[, keep], meta_s[keep, ], "HSCT_responder")

# Fisher test
T52 <- prevalence(res$presence, res$absence)
T52[is.na(T52)] <- 1
T52 <- p.adjust(T52, "BH")

any(T52 < 0.05)

# Test the prevalence between non controls at times T0, T26, T52 ####
removeControls <- meta_s$Controls  == "NO"
removeTimes <- meta_s$Time %in% c("T0", "T26", "T52")
keep <- removeControls & removeTimes
prevalence_d <- sweep(otus_table_s[, keep], 2, colSums(otus_table_s[, keep]), `/`) > 0.005
subSets <- allComb(meta_s[keep, ], c("Time", "HSCT_responder"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]
presence <- prevalence_d %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                  byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
IBD_Time <- prevalence(presence, absence)

# Remove and adjust
IBD_Time[is.na(IBD_Time)] <- 1
IBD_Time <- p.adjust(IBD_Time, "BH")

if (any(IBD_Time < 0.05)) {
  message("Plots!!")
}
# The presence of microorganisms thorough time is not different in the 
# overall patients

# Test if the responders and the non responders behave differently along time
# If the presence of microorganisms thorough time is different in responders 
# than in non-responders ####
removeControls <- meta_s$Controls  == "NO"
keepResponders <- meta_s$HSCT_responder[removeControls] == "YES"
keepNonResponders <- meta_s$HSCT_responder[removeControls] == "NO"
wocontrols <- otus_table_s[, removeControls]

subSets <- allComb(meta_s[removeControls, ][keepResponders, ], c("Time"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]

prevalence_d <- sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005

presence <- prevalence[, keepResponders] %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                        byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
responders <- prevalence(presence, absence)


subSets <- allComb(meta_s[removeControls, ][keepNonResponders, ], c("Time"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]

presence <- prevalence_d[, keepNonResponders] %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                        byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
nonResponders <- prevalence(presence, absence)


# Bootstrapp while controlling by Time ####
# bootstrapping with 1000 replications
# results <- boot(data = t(sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005),
#                 statistic = ratio,
#                 R = 1000,
#                 meta = meta_s[removeControls, ], columns = "Time",
#                 parallel = "multicore", ncpus = 4)
# 
# # view results
# results
# plot(results)
# 
# # get 95% confidence interval
# boot.ci(results, type="bca")

dev.off()