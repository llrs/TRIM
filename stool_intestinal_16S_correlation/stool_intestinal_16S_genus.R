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
eqOTUS <- read.csv("equivalent_genus.csv", stringsAsFactors = FALSE)
setwd(cd)

# Summarize to genus
library("metagenomeSeq")


# Create the objects
MR_i <- newMRexperiment(t(otus_i), 
                      phenoData = AnnotatedDataFrame(meta),
                      featureData = AnnotatedDataFrame(as.data.frame(tax_i)))
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")

MR_s <- newMRexperiment(t(otus_s), 
                        phenoData = AnnotatedDataFrame(meta),
                        featureData = AnnotatedDataFrame(as.data.frame(tax_s)))
genus_s <- aggTax(MR_s, lvl = "Genus", out = "matrix")


# Identify the species
ta <- tax_i[eqOTUS$intestinal, "Genus", drop = FALSE]
ta$Genus <- factor(ta$Genus)
tas <- ta$Genus
names(tas) <- rownames(ta)

## OTUs of the species in both datasets
u_i <- unique(eqOTUS$intestinal)
u_s <- unique(eqOTUS$stools)

# Known which OTUs is which specie
eqOTUS$Genus <- tas[eqOTUS$intestinal]
for (s in eqOTUS$stools) {
  sp <- eqOTUS$Genus[eqOTUS$stools == s]
  eqOTUS$Genus[eqOTUS$stools == s] <- sp[!is.na(sp)]
}

# Colors for the different OTUs of the same Species
colIntestinal <- split(eqOTUS$Genus, eqOTUS$intestinal)
colIntestinal <- sapply(colIntestinal, function(x)x[1])

colStools <- split(eqOTUS$Genus, eqOTUS$stools)
colStools <- sapply(colStools, function(x)x[1])

# Set the same color for the microorganism for different OTUs of the different
# datasets
palette(rainbow(length(unique(tas))))
color_species <- palette()
color_species
levels(colIntestinal) <- color_species
levels(colStools) <- color_species

# Set the order of the levels in Time
meta$Time <- factor(as.character(meta$Time), 
                    levels = c("T0", "T26", "T52", "T106", "TM36", "TM48"))

corOTUS <- function(x){
  if (sum(x, na.rm = TRUE) == 0) {
    return(NA)
  }
  o <- cor(otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs", 
           method = "spearman")
  colNA <- apply(is.na(o), 2, all)
  rowNA <- apply(is.na(o), 1, all)
  o2 <- o[!rowNA, !colNA]
  
  colSide <- colnames(o2)
  rowSide <- rownames(o2)
  colorC <- order(colStools[colSide])
  colorR <- order(colIntestinal[rowSide], decreasing = TRUE)
  heatmap(o2[colorR, colorC], Rowv = NA, Colv = NA, scale = "none", 
                ColSideColors = as.character(colStools[colSide][colorC]),
                RowSideColors = as.character(colIntestinal[rowSide][colorR]), 
                xlab = "Stools OTUs", ylab = "Biopsies OTUs")
  o
}

cors2OTUs <- function(y){
  
  # Find all the OTUs in the intestinal
  o_i <- eqOTUS$intestinal[eqOTUS$stools == y]
  # Find all the otus of stools for the same org. 
  o_s <- lapply(o_i, function(x){eqOTUS$stools[eqOTUS$intestinal == x]})
  o_s <- unique(unlist(o_s))
  
  if (!y %in% o_s) {
    stop(y, " should be on the list of otus ", paste(o_i, collapse = ", "))
  }
  
  if (is.list(subCors) && length(subCors) >= 2 ) {
    test <- sapply(seq_along(subCors), function(xy) {
      if (is.na(subCors[[xy]])) {
        return(NA)
      }
      mOTUs <- mean(subCors[[xy]][o_i, o_s], na.rm = TRUE)
      
      names(mOTUs) <- names(subCors)[xy]
      mOTUs
    })
  } else {
    test <- mean(subCors[o_i, o_s], na.rm = TRUE)
  }
  
  test
}

pdf(paste0("Figures/", today, "_correlations_genus.pdf"))

h <- cor(t(genus_i), t(genus_s))

heatmap(h, scale = "none", main = "Correlations between Genus in all samples")

cors <- corOTUS(rep(TRUE, nrow(otus_i)))

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
  
  ta <- tax_s[y, "Genus"]
  mCors <- mean(cors[o_i, o_s], na.rm = TRUE)
  c("All" = mCors, "ta" = ta)
})
corsOrg <- as.data.frame(t(corsOrg))
corsOrg$All <- as.numeric(levels(corsOrg$All))[corsOrg$All] 
ggplot(corsOrg, aes("All", y = All)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") +
  ylab("Correlation") +
  ylim(-1, 1)


### Compare the equivalent otus in different settings
area <- allComb(meta, c("CD_Aftected_area"))
subCors <- sapply(as.data.frame(area), corOTUS, simplify = FALSE)

# Do the mean of those correlations that correspond to the same microorganism
# before doing the boxplots
cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), cors2OTUs)

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE) +
  ylim(-1, 1)

corsOrg <- cbind(corsOrg, t(cors2Org))
ggplot(melt(corsOrg), aes(variable, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Samples") + 
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(corsOrg, aes("ratio", y = COLON/ILEUM)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(aes(col = ta), position = position_jitter()) +
  xlab("Samples") +
  guides(col = FALSE)

### Compare the equivalent otus in different settings
time <- allComb(meta, "Time")
subCors <- sapply(as.data.frame(time), corOTUS, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), cors2OTUs)
corsOrg <- cbind(corsOrg, t(cors2Org))

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Time") + 
  ylab("Correlation rho") + 
  guides(col = FALSE) +
  ylim(-1, 1)

responders <- allComb(meta, "HSCT_responder")
subCors <- sapply(as.data.frame(responders), corOTUS, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), cors2OTUs)
corsOrg <- cbind(corsOrg, t(cors2Org))

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Responder") + 
  ylab("Correlation rho") +
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(corsOrg, aes("Ratio responders", y = YES/NO)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Responder") + 
  guides(col = FALSE)

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time", "CD_Aftected_area"))
subCors <- sapply(as.data.frame(time_area), corOTUS, simplify = FALSE)

cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), cors2OTUs)
corsOrg <- cbind(corsOrg, t(cors2Org))

ggplot(melt(t(cors2Org)), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Division by time and afected area") + 
  ylab("Correlation rho") +
  guides(col = FALSE) +
  ylim(-1, 1)

### Compare the equivalent otus in different settings
time_responder <- allComb(meta, c("Time", "HSCT_responder"))
subCors <- sapply(as.data.frame(time_responder), corOTUS, simplify = FALSE)
cors2Org <- sapply(rownames(unique(tax_s[eqOTUS$stools,])), cors2OTUs)
cors2Org <- t(cors2Org)
keep <- apply(cors2Org, 2, function(x){all(is.na(x))})
cors2Org <- cors2Org[, !keep]

ggplot(melt(cors2Org), aes(Var2, y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  # geom_point(position = position_jitter()) +
  xlab("Division by time and responders") + 
  ylab("Correlation rho") + 
  guides(col = FALSE) +
  ylim(-1, 1)

corsOrg <- cbind(corsOrg, cors2Org)


## 22 more common organisms when we consider "Time", "Endoscopic_Activity", 
## "Treatment", "CD_Aftected_area", "Involved_Healthy", "Active_area", "ID"

# 
comorg <- read.csv("../stool_intestinal_metadb/important_common_microrg.csv", 
                   stringsAsFactors = FALSE)

ta_common <- unique(comorg[comorg$Genus != "", "Genus"])
corsOrg$Important <- corsOrg$ta %in% ta_common
ggplot(melt(corsOrg[, c("Important", "All")]), aes(Important, y = value)) + 
  # geom_violin() +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 1)) +
  xlab("Important species") + 
  ylab("Correlation rho") +
  ggtitle("Correlation in all samples") +
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(melt(corsOrg[, c("Important", "COLON")]), aes("COLON", y = value)) + 
  # geom_violin() +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 1)) +
  xlab("Important species") + 
  ylab("Correlation rho") +
  ggtitle("Correlation in colon samples with their stools") +
  facet_wrap(~Important) +
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(melt(corsOrg[, c("Important", "ILEUM")]), aes("ILEUM", y = value)) + 
  # geom_violin() +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 1)) +
  xlab("Important species") + 
  ylab("Correlation rho") +
  ggtitle("Correlation in ileum samples with their stools") +
  facet_wrap(~Important) +
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(melt(corsOrg[, c("Important", "YES")]), aes("Responder", y = value)) + 
  # geom_violin() +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 1)) +
  xlab("Important species") + 
  ylab("Correlation rho") +
  ggtitle("Correlation in responders") +
  facet_wrap(~Important) +
  guides(col = FALSE) +
  ylim(-1, 1)

ggplot(melt(corsOrg[, c("Important", "NO")]), aes("Non responders", y = value)) + 
  # geom_violin() + 
  geom_boxplot() +
  geom_point(position = position_dodge(width = 1)) +
  xlab("Important species") + 
  ylab("Correlation rho") +
  ggtitle("Correlation in Non responders") +
  facet_wrap(~Important) +
  guides(col = FALSE) +
  ylim(-1, 1)

dev.off() 