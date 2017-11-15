cd <- setwd("..")

# Load the helper file
source("helper_functions.R")
# Read files
otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)
tax_i <- read.csv(file = "intestinal_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)
tax_s <- read.csv(file = "stools_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)
eqOTUS <- read.csv("equivalent_otus.csv", stringsAsFactors = FALSE)
setwd(cd)

keep <- !grepl("28_T52_T_DM_CH", meta$Sample_Code) # Remove outlier See PCA

meta <- meta[keep, ]
otus_s <- otus_s[keep, ]
otus_i <- otus_i[keep, ]

# Identify the species
ta <- unique(tax_i[eqOTUS$intestinal, c("Genus", "Species")])
ta <- apply(ta, 1, paste, collapse = " ")

## OTUs of the species in both datasets
u_i <- unique(eqOTUS$intestinal)
u_s <- unique(eqOTUS$stools)

findCorOTUs <- function(x){
  if (sum(x) == 0) {
    return(NA)
  } else {
  cor(otus_i[x, u_i], otus_s[x, u_s], 
      use = "pairwise.complete.obs", method = "spearman")
  }
}
pdf(paste0("Figures/", today, "_correlations.pdf"))
cors <- findCorOTUs(rep(TRUE, nrow(otus_i)))

# Do the mean of those correlations that correspond to the same microorganism
# before doing the boxplots
corsOrg <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  mCors <- mean(cors[o_i, o_s], na.rm = TRUE)
  c("All" = mCors, "ta" = ta)
})
corsOrg <- as.data.frame(t(corsOrg))
corsOrg$All <- as.numeric(levels(corsOrg$All))[corsOrg$All] 
ggplot(corsOrg, aes("All", y = All)) + 
  geom_violin() + 
  # geom_boxplot() +
  geom_point(position = position_jitter()) +
  xlab("Samples")
  

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("CD_Aftected_area"))
subCors <- sapply(as.data.frame(time_area), findCorOTUs, simplify = FALSE)

# Do the mean of those correlations that correspond to the same microorganism
# before doing the boxplots
cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    
    names(mOTUs) <- names(subCors)[xy]
    mOTUs
  })
  test
})

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

corsOrg <- cbind(corsOrg, t(cors2Org))
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

ggplot(corsOrg, aes("ratio", y = COLON/ILEUM)) + 
  geom_violin() + 
  # geom_boxplot() +
  geom_point(aes(col = ta), position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

### Compare the equivalent otus in different settings
time <- allComb(meta, "Time")
subCors <- sapply(as.data.frame(time), function(x){
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    
    names(mOTUs) <- names(subCors)[xy]
    mOTUs
  })
  test
})
corsOrg <- cbind(corsOrg, t(cors2Org))
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)
time <- allComb(meta, "HSCT_responder")
subCors <- sapply(as.data.frame(time), function(x){
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    
    names(mOTUs) <- names(subCors)[xy]
    mOTUs
  })
  test
})
corsOrg <- cbind(corsOrg, t(cors2Org))
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time", "CD_Aftected_area"))
subCors <- sapply(as.data.frame(time_area), function(x){
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    
    names(mOTUs) <- names(subCors)[xy]
    mOTUs
  })
  test
})
corsOrg <- cbind(corsOrg, t(cors2Org))
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

### Compare the equivalent otus in different settings
time_responder <- allComb(meta, c("Time", "HSCT_responder"))
subCors <- sapply(as.data.frame(time_responder), function(x){
  if (sum(x) == 0) {
    return(NA)
  }
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)
cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    if (is.null(dim(subCors[[xy]]))) {
      return(NA)
    }
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    
    names(mOTUs) <- names(subCors)[xy]
    mOTUs
  })
  test
})
cors2Org <- t(cors2Org)
keep <- apply(cors2Org, 2, function(x){all(is.na(x))})
cors2Org <- cors2Org[, !keep]

ggplot(melt(cors2Org), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

corsOrg <- cbind(corsOrg, cors2Org)
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE)

dev.off() 
save.image(file = paste0(today, "_stool_intestinal.RData"))