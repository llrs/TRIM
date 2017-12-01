library("metagenomeSeq")
library("phyloseq")

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

pdf(paste0("Figures/", today, "_plots.pdf"))

# Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x){length(unique(x)) != 1})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)

# Create the file
MR <- newMRexperiment(otus_table_s[, ordSamples], 
                      phenoData = AnnotatedDataFrame(meta_s[ordSamples, ]),
                      featureData = AnnotatedDataFrame(as.data.frame(otus_tax_s)))
filterData(MR, present = 10, depth = 1000)
p <- cumNormStatFast(MR)
MR <- cumNorm(MR, p = p)
assayData(MR)$relative <- assayData(MR)$counts/rowSums(assayData(MR)$counts)
assayData(MR)$prevalence <- as.matrix(assayData(MR)$counts != 0)

rareFeatures <- which(rowSums(MRcounts(MR) > 0) < 10)
MRtrim <-  MR[-rareFeatures, ]
MRp <- cumNormStat(MRtrim, pFlag = TRUE, main = "Trimmed lung data")
MRtrim <- cumNorm(MRtrim, p = MRp)
assayData(MRtrim)



normFactor <- normFactors(MRtrim)
assayData(MRtrim)$relative <- assayData(MRtrim)$counts/rowSums(assayData(MRtrim)$counts)
assayData(MRtrim)$prevalence <- as.matrix(assayData(MRtrim)$counts != 0)

normFactor <- log2(normFactor/median(normFactor) + 1)
settings <- zigControl(maxit = 10, verbose = TRUE)
mod <- model.matrix(~ 0 + Endoscopic_Activity, data = pData(MRtrim))
mod <- cbind(mod, ID = pData(MRtrim)$ID, Time = as.numeric(as.factor(pData(MRtrim)$Time)))
fit <- fitZig(obj = MRtrim, mod = mod, useCSSoffset = FALSE,
              control = settings)

eb <- fit$eb
head(cbind(eb$F, eb$F.p.value))
# From now on like with limma
# To consider if we want to test for specific phylum or taxa 
?metagenomeSeq::aggregateByTaxonomy

dev.off()