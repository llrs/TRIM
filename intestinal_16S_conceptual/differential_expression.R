library("metagenomeSeq")
library("phyloseq")

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
ordSamples <- intersect(rownames(meta_i), colnames(otus_table_i))

# Create the object
MR <- newMRexperiment(otus_table_i[, ordSamples], 
                      phenoData = AnnotatedDataFrame(meta_i[ordSamples, ]),
                      featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i)))
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