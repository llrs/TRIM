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


## species that are in both datasets
u_i <- unique(eqOTUS$intestinal)
u_s <- unique(eqOTUS$stools)

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("CD_Aftected_area"))
subCors <- sapply(as.data.frame(time_area), function(x){
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)

# Tranform the data for the plot
cors <- sapply(subCors, function(x){
  x <- x[upper.tri(x)]
  x[!is.na(x)]
})

cors <- cors[colSums(time_area) >= 4]
corEqOtus <- melt(cors)

# Set the factors of Time in the order we want
ggplot(corEqOtus) +
  geom_boxplot(aes(L1, value)) + 
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools") + 
  xlab("CD afected area")

# Do the mean of those correlations that correspond to the same microorganism
# before doing the boxplots
a <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), function(y){
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  o_s <- lapply(o_i, function(x){
    eqOTUS$stools[eqOTUS$intestinal == x]
  })
  o_s <- unique(unlist(o_s))
  o_i <- lapply(o_s, function(x){
    eqOTUS$intestinal[eqOTUS$stools == x]
  })
  o_i <- unique(unlist(o_i))
  
  d <- matrix(FALSE, ncol = ncol(subCors[[1]]), nrow = nrow(subCors[[1]]), 
              dimnames = dimnames(subCors[[1]]))
  d[o_i, o_s] <- TRUE
  ta <- paste(tax_s[y, c("Genus", "Species")], collapse = " ")
  test <- sapply(seq_along(subCors), function(xy) {
    # heatmap(t(subCors[[xy]]), scale = "none", 
    #         ylab = "OTUs stools", xlab = "OTUs intestinal",
    #         Rowv = NA, Colv = NA,  
    #         main = paste0("Correlation of ", ta, " in ", names(subCors)[xy]), 
    #         add.expr = { makeRects(d,"blue")})
    mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
    list("Species" = ta, "mean_cor" = mOTUs, "cond" = names(subCors)[xy], 
         "n" = sum(time_area[, xy], na.rm = TRUE))
  })
  t(test)
})


### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time", "CD_Aftected_area"))
subCors <- sapply(as.data.frame(time_area), function(x){
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)


# Tranform the data for the plot
cors <- sapply(subCors, function(x){
  x <- x[upper.tri(x)]
  x[!is.na(x)]
})

cors <- cors[colSums(time_area) >= 4]
corEqOtus <- melt(cors)

name <- strsplit(corEqOtus$L1, "_|_")
corEqOtus$Time <- sapply(name, function(x){x[1]})
corEqOtus$Site <- sapply(name, function(x){x[3]})
# Set the factors of Time in the order we want
corEqOtus$Time <- factor(corEqOtus$Time, 
                         levels = c("T0", "TM36", "TM48", "T26", "T52", "T106"))

ggplot(corEqOtus) +
  geom_boxplot(aes(Site, value)) + 
  facet_grid(~Time) +
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools")

### Compare the equivalent otus in different settings
time_responder <- allComb(meta, c("Time", "HSCT_responder"))
subCors <- sapply(as.data.frame(time_responder), function(x){
  if (sum(x) == 0) {
    return(NA)
  }
  cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
      method = "spearman")
}, simplify = FALSE)


# Tranform the data for the plot
cors <- sapply(subCors, function(x){
  x <- x[upper.tri(x)]
  x[!is.na(x)]
})

cors <- cors[colSums(time_responder) >= 4]
corEqOtus <- melt(cors)

name <- strsplit(corEqOtus$L1, "_|_")
corEqOtus$Time <- sapply(name, function(x){x[1]})
corEqOtus$Responder <- sapply(name, function(x){x[3]})
# Set the factors of Time in the order we want
corEqOtus$Time <- factor(corEqOtus$Time, 
                         levels = c("T0", "TM36", "TM48", "T26", "T52", "T106"))

ggplot(corEqOtus) +
  geom_boxplot(aes(Responder, value)) + 
  facet_grid(~Time) +
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools")