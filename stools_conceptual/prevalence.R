library("boot")

wd <- setwd("..")
source("helper_functions.R")
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

pdf(paste0("Figures/", today, "_ratios_plots.pdf"))

 # Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x){length(unique(x)) != 1})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)

# Test the prevalence between controls and non controls ####
prevalence <- sweep(otus_table_s, 2, colSums(otus_table_s), `/`) > 0.005
subSets <- allComb(meta_s, c("Controls"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]
presence <- prevalence %*% subSets
absence <- matrix(totalSamples, nrow = nrow(presence), 
                   byrow = TRUE, ncol = ncol(presence)) - presence
# Fisher test
d <- sapply(rownames(presence), function(i){
  m <- rbind(P = presence[i, ], 
             A = absence[i, ])
  f <- fisher.test(m, workspace = 2e8)
  # ch <- chisq.test(m)
  # c(f$p.value, ch$p.value)
  f$p.value
})
d[is.na(d)] <- 1
d <- p.adjust(d, "BH")
if (any(d < 0.05)) {
  ta <- unique(otus_tax_s[d < 0.05, c("Genus", "Species")])
  species <- unique(ta[, 2])
  species[!is.na(species)]
  
  # Plot the 
  genus <- ta[, 1]
  genus <- t(t(table(genus)))
  genus <- cbind(genus, rownames(genus))
  colnames(genus) <- c("OTUs", "Genus")
  genus <- as.data.frame(genus)
  genus[, "OTUs"] <- as.numeric(genus[, "OTUs"])
  genus <- genus[order(genus$OTUs, decreasing = TRUE), ]
  genus$Genus <- as.character(genus$Genus)
  genus$Genus <- factor(genus$Genus, levels = genus$Genus)
  
  ggplot(genus) + 
    geom_bar(aes(Genus, OTUs), stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Genus related to having IBD")
}
# summary:  Not significative, IBD patients are as  likely to have one 
# microorganism as the controls.

# Test the prevalence between non controls at times T0, T26, T52 ####
removeControls <- meta_s$Controls  == "NO"
removeTimes <- !meta_s$Time %in% c("T0", "T26", "T52")
keep <- removeControls & removeTimes
prevalence <- sweep(otus_table_s[, keep], 2, colSums(otus_table_s[, keep]), `/`) > 0.005
subSets <- allComb(meta_s[keep, ], c("Time", "HSCT_responder"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]
presence <- prevalence %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                  byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
d <- sapply(rownames(presence), function(i) {
  m <- rbind(P = presence[i, ], 
             A = absence[i, ])
  f <- fisher.test(m, workspace = 2e8)
  # ch <- chisq.test(m)
  # c(f$p.value, ch$p.value)
  f$p.value
})
d[is.na(d)] <- 1
d <- p.adjust(d, "BH")
if (any(d < 0.05)) {
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

prevalence <- sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005

presence <- prevalence[, keepResponders] %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                        byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
responders <- sapply(rownames(presence), function(i) {
  m <- rbind(P = presence[i, ], 
             A = absence[i, ])
  f <- fisher.test(m, workspace = 2e8)
  # ch <- chisq.test(m)
  # c(f$p.value, ch$p.value)
  f$p.value
})

subSets <- allComb(meta_s[removeControls, ][keepNonResponders, ], c("Time"))
totalSamples <- colSums(subSets)
subSets <- subSets[, totalSamples >= 1]
totalSamples <- totalSamples[totalSamples >= 1]

presence <- prevalence[, keepNonResponders] %*% subSets
totalSamplesm <- matrix(totalSamples, nrow = nrow(presence), 
                        byrow = TRUE, ncol = ncol(presence)) 
absence <- totalSamplesm - presence

# Fisher test
nonResponders <- sapply(rownames(presence), function(i) {
  m <- rbind(P = presence[i, ], 
             A = absence[i, ])
  f <- fisher.test(m, workspace = 2e8)
  # ch <- chisq.test(m)
  # c(f$p.value, ch$p.value)
  f$p.value
})


# Bootstrapp while controlling by Time ####
ratio <- function(columns, data, indices, meta) {
  a <- t(data[indices, ]) # allows boot to select sample
  bindices <- !seq_len(nrow(data)) %in% indices
  b <- t(data[bindices, ])
  Ameta <- meta[indices, ]
  Bmeta <- meta[bindices, ]
  
  AsubSets <- allComb(Ameta, columns)
  if (!is.matrix(AsubSets)) {
    return(NA)
  }
  AtotalSamples <- colSums(AsubSets)
  AsubSets <- AsubSets[, AtotalSamples >= 1]
  AtotalSamples <- AtotalSamples[AtotalSamples >= 1]
  
  Apresence <- a %*% AsubSets
  AtotalSamplesm <- matrix(AtotalSamples, nrow = nrow(Apresence), 
                          byrow = TRUE, ncol = ncol(Apresence)) 
  Aabsence <- AtotalSamplesm - Apresence
  
  BsubSets <- allComb(Bmeta, columns)
  if (!is.matrix(BsubSets)) {
    return(NA)
  }
  BtotalSamples <- colSums(BsubSets)
  BsubSets <- BsubSets[, BtotalSamples >= 1]
  BtotalSamples <- BtotalSamples[BtotalSamples >= 1]
  
  Bpresence <- b %*% BsubSets
  BtotalSamplesm <- matrix(BtotalSamples, nrow = nrow(Bpresence), 
                          byrow = TRUE, ncol = ncol(Bpresence)) 
  Babsence <- BtotalSamplesm - Bpresence
  
  
  
  # Fisher test and ratio calculation
  sapply(rownames(presence), function(i) {
    Am <- rbind(P = Apresence[i, ], 
               A = Aabsence[i, ])
    Af <- fisher.test(Am, workspace = 2e8)
    
    Bm <- rbind(P = Bpresence[i, ], 
                A = Babsence[i, ])
    Bf <- fisher.test(Bm, workspace = 2e8)
    
    # ch <- chisq.test(m)
    # c(f$p.value, ch$p.value)
    Af$p.value/Bf$p.value
  })
}
# 
# # bootstrapping with 1000 replications
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