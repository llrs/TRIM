intestinal <- "intestinal_16S"
stool <- "stools_16S"

library("ggplot2")
source("helper_functions.R")

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Read the stools OTUs
otus_table_s <- read.delim(file.path(stool, "OTUs-Table-refined-stools.tab"), 
                           stringsAsFactors = FALSE, row.names = 1,
                           check.names = FALSE)
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

theme_set(theme_bw())
# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(file_meta_s, check.names = FALSE, row.names = 1, 
                     stringsAsFactors = FALSE)
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)
# Make PCAs
keep_i <- !grepl("28_T52_T_DM_CH", meta_i$Sample_Code) # Remove outlier
pca_i <- prcomp(t(otus_table_i[, keep_i]), scale. = TRUE)
pca_s <- prcomp(t(otus_table_s), scale. = TRUE)

pca_s_x <- as.data.frame(pca_s$x)
pca_i_x <- as.data.frame(pca_i$x)

ggplot(pca_s_x) +
  geom_text(aes(PC1, PC2, col = meta_s$Patient_ID, 
                 label = meta_s$Time)) +
  ggtitle("PCA stools") + 
  guides(col = guide_legend(title="Patient ID"))

label <- strsplit(as.character(meta_i$Sample_Code)[keep_i], split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[2], x[5], sep = "_")}
  else if (length(x) != 5) {
    paste(x[1], x[4], sep = "_")
    }
  })


ggplot(pca_i_x) +
  geom_text(aes(PC1, PC2, col = meta_i$Patient_ID[keep_i], label = labels)) +
  ggtitle("PCA mucosa") + 
  guides(col = guide_legend(title="Patient ID"))
