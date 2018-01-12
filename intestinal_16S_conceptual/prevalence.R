library("boot")

wd <- setwd("..")
source("helper_functions.R")
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


# Test the prevalence between controls and non controls ####
prevalence <- sweep(otus_table_i, 2, colSums(otus_table_i), `/`) > 0.005
subSets <- allComb(meta_i, c("IBD"))
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
removeControls <- meta_i$IBD  == "CD"
removeTimes <- !meta_i$Time %in% c("T0", "T26", "T52")
keep <- removeControls & removeTimes
prevalence <- sweep(otus_table_i[, keep], 2, colSums(otus_table_i[, keep]), `/`) > 0.005
subSets <- allComb(meta_i[keep, ], c("Time", "HSCT_responder"))
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
removeControls <- meta_i$IBD  == "CD"
keepResponders <- meta_i$HSCT_responder[removeControls] == "YES"
keepNonResponders <- meta_i$HSCT_responder[removeControls] == "NO"
wocontrols <- otus_table_i[, removeControls]

subSets <- allComb(meta_i[removeControls, ][keepResponders, ], c("Time"))
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

subSets <- allComb(meta_i[removeControls, ][keepNonResponders, ], c("Time"))
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

pdf(paste0("Figures/", today, "_ratios_plots.pdf"))
hist(log10(nonResponders/responders))

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

# bootstrapping with 1000 replications
results <- boot(data = t(sweep(wocontrols, 2, colSums(wocontrols), `/`) > 0.005),
                statistic = ratio,
                R = 1000,
                meta = meta_i[removeControls, ], columns = "Time",
                parallel = "multicore", ncpus = 4)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")

dev.off()