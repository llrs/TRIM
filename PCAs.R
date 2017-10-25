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
pca_s <- prcomp(t(otus_table_s), scale. = TRUE)

pca_s_x <- as.data.frame(pca_s$x)

ggplot(pca_s_x) +
  geom_text(aes(PC1, PC2, col = meta_s$Patient_ID, label = meta_s$Time)) +
  ggtitle("PCA stools") + 
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5))


names(colors) <- unique(meta_i$Patient_ID)

time <- "T0"

keepSub <- meta_i$Time == time & !grepl("28_T52_T_DM_CH", meta_i$Sample_Code)
subOtu <- otus_table_i[, keepSub]

pca_i <- prcomp(t(subOtu[rowSums(subOtu) != 0, ]), 
                scale. = TRUE)
subMeta_i <- meta_i[keepSub, ]
pca_i_x <- as.data.frame(pca_i$x)

label <- strsplit(as.character(subMeta_i$Sample_Code), 
                  split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[1], x[5], sep = "_")
    # x[5]
    }
  
  else if (length(x) != 5) {
    paste(x[1], x[4], sep = "_")
    # x[4]
    }
  })

ggplot(pca_i_x) +
  geom_text(aes(PC1, PC2, col = meta_i$Patient_ID[keepSub], 
                label = labels)) + 
  scale_color_manual(values = colors) +
  ggtitle(paste0("PCA mucosa at time (", time, ")")) + 
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5))

label <- strsplit(as.character(meta_i$Sample_Code)[keep_i], 
                  split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[1], x[2], x[5], sep = "_")
    # x[5]
  }
  
  else if (length(x) != 5) {
    paste(x[1], x[4], sep = "_")
    # x[4]
  }
})
ggplot(pca_i_x) +
  geom_text(aes(PC1, PC2, col = meta_i$Patient_ID[keep_i], 
                label = labels)) + 
  scale_color_manual(values = colors) +
  ggtitle("PCA mucosa ") + 
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(range(pca_i_x[, 2])) + 
  xlim(range(pca_i_x[, 1]))



# Count with a cutofff of 00.5 in the relative abundance!!!
meta_ii <- meta_i[!meta_i$Patient_ID %in% c("40", "41", "38") 
                  & grepl("INVOLVED", meta_i$Involved_Healthy) 
                  & meta_i$Time %in% c("T0", "T26", "T52") |
                    grepl("^C", meta_i$Patient_ID), ]
lacto <- lacto[, !meta_i$Patient_ID %in% c("40", "41", "38") 
              & grepl("INVOLVED", meta_i$Involved_Healthy) 
              & meta_i$Time %in% c("T0", "T26", "T52") |
                grepl("^C", meta_i$Patient_ID)]
sum(colSums(lacto[, meta_ii$Time == "C" & 
                    # meta_ii$HSCT_responder == "YES" &
                    # meta_ii$ == "INVOLVED" & 
                    grepl("COLON", meta_ii$CD_Aftected_area)]) != 0)

# sum(colSums(lacto[, meta_i$Time == "T52" & 
                    meta_i$HSCT_responder == "YES" & 
                    meta_i$Involved_Healthy == "INVOLVED" & 
                    grepl("COLON", meta_i$CD_Aftected_area)]) != 0)