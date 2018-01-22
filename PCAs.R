intestinal <- "intestinal_16S"
stool <- "stools_16S"
rna <- "intestinal_RNAseq"

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

# Load the RNAseq
expr <- read.delim(file.path(rna, "table.counts.results"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(file_meta_s, check.names = FALSE, row.names = 1, 
                     stringsAsFactors = FALSE)
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)

file_meta_r <- file.path(rna, "111217_metadata.csv")
meta_r <- read.table(file_meta_r, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = ";", 
                     na.strings = c("NA", ""))
colnames(meta_r) <- meta_r[1, ]
meta_r <- meta_r[-1, ]

# Clean the metadata
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i <- meta_i[, apply(meta_i, 2, function(x){length(unique(x)) != 1})]
meta_i$Active_area[meta_i$Active_area == ""] <- NA
meta_i$ID <- meta_i$Patient_ID
meta_i$ID[meta_i$Patient_ID %in% c("15", "23")] <- "15/23"
meta_i$ID[meta_i$Patient_ID %in% c("33", "36")] <- "33/36"
meta_i$ID[meta_i$Patient_ID %in% c("29", "35")] <- "29/35"
meta_i$ID <- as.factor(meta_i$ID)
# There is a mislabeling on those tubes, we don't know which is which
meta_i$CD_Aftected_area[meta_i$Sample_Code == "22_T52_T_DM_III"] <- NA


# Clean the metadata
meta_s <- meta_s[, apply(meta_s, 2, function(x){length(unique(x)) != 1})]
meta_s$ID <- meta_s$Patient_ID
meta_s$ID[meta_s$Patient_ID %in% c("15", "23")] <- "15/23"
meta_s$ID[meta_s$Patient_ID %in% c("33", "36")] <- "33/36"
meta_s$ID[meta_s$Patient_ID %in% c("29", "35")] <- "29/35"
meta_s$ID <- as.factor(meta_s$ID)

# Clean the metadata
meta_r$Active_area[meta_r$Active_area == ""] <- NA
meta_r <- meta_r[, apply(meta_r, 2, function(x){length(unique(x)) != 1})]
meta_r$Active_area[meta_r$Active_area == ""] <- NA
meta_r$ID <- meta_r$Patient_ID
meta_r$ID[meta_r$Patient_ID %in% c("15", "23")] <- "15/23"
meta_r$ID[meta_r$Patient_ID %in% c("33", "36")] <- "33/36"
meta_r$ID[meta_r$Patient_ID %in% c("29", "35")] <- "29/35"
meta_r$ID <- as.factor(meta_r$ID)

# Make PCAs
pca_s <- prcomp(t(otus_table_s), scale. = TRUE)
pca_s_x <- as.data.frame(pca_s$x)
pca_s_var <- round(summary(pca_s)$importance[2, ]*100, digits = 2)


# Define colors
colors_i <- colors
colors_s <- colors
colors_ir <- colors
names(colors_i) <- unique(meta_i$ID)
names(colors_s) <- unique(meta_s$ID)
names(colors_ir) <- unique(meta_r$ID)


pcas <- cbind(pca_s_x, meta_s)

pdf(paste0("Figures/", today, "_PCA.pdf"))


ggplot(pcas) +
  geom_text(aes(PC1, PC2, col = ID, label = Time)) +
  ggtitle("PCA stools") + 
  guides(col = guide_legend(title = "Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors_s) +
  xlab(paste("PC1", pca_s_var[1], "%")) +
  ylab(paste("PC2", pca_s_var[2], "%"))

ggplot(pcas) +
  geom_text(aes(PC1, PC2, col = HSCT_responder, 
                label = paste(Time, ID, sep = "_"))) +
  ggtitle("PCA stools") + 
  guides(col = guide_legend(title = "Responders")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_s_var[1], "%")) +
  ylab(paste("PC2", pca_s_var[2], "%"))

# PCA intestinals

pca_i <- prcomp(t(otus_table_i), scale. = TRUE)
pca_i_x <- as.data.frame(pca_i$x)
pca_i_var <- round(summary(pca_i)$importance[2, ]*100, digits = 2)

label <- strsplit(as.character(meta_i$Sample_Code), split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[5], sep = "_")
    # x[5]
  }
  
  else if (length(x) != 5) {
    paste(x[4], sep = "_")
    # x[4]
  }
})
meta_i <- cbind(meta_i, labels)

pcai <- cbind(pca_i_x, meta_i)

ggplot(pcai) +
  geom_text(aes(PC1, PC2, col = ID, label = labels)) + 
  guides(col = guide_legend(title = "Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors_i) +
  xlab(paste("PC1", pca_i_var[1], "%")) +
  ylab(paste("PC2", pca_i_var[2], "%")) + 
  ggtitle("PCA biopsies")

# Remove the patient We can see that there is out of the line
keep <- !grepl("28_T52_T_DM_CH", meta_i$Sample_Code)
meta_i <- meta_i[keep, ]
otus_table_i <- otus_table_i[, keep]

pca_i <- prcomp(t(otus_table_i), scale. = TRUE)
pca_i_x <- as.data.frame(pca_i$x)
pca_i_var <- round(summary(pca_i)$importance[2, ]*100, digits = 2)

pcai <- cbind(pca_i_x, meta_i)

ggplot(pcai) +
  geom_text(aes(PC1, PC2, col = ID, label = labels)) + 
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors_i) +
  xlab(paste("PC1", pca_i_var[1], "%")) +
  ylab(paste("PC2", pca_i_var[2], "%")) + 
  ggtitle("PCA biopsies")

ggplot(pcai) +
  geom_text(aes(PC1, PC2, col = HSCT_responder, 
                label = paste(Time, ID, sep = "_"))) +
  ggtitle("PCA biopsies") + 
  guides(col = guide_legend(title = "Responders")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_i_var[1], "%")) +
  ylab(paste("PC2", pca_i_var[2], "%"))


# PCA intestinal RNAseq with Barcelona
# Filter expression
expr <- norm_RNAseq(expr)

pca_ir <- prcomp(t(expr), scale. = TRUE)
pca_ir_x <- as.data.frame(pca_ir$x)
pca_ir_var <- round(summary(pca_ir)$importance[2, ]*100, digits = 2)

pca_ir_x <- pca_ir_x[rownames(pca_ir_x) %in% meta_r$`Sample Name_RNA`, ]
meta_r_ord <- meta_r[meta_r$`Sample Name_RNA` %in% rownames(pca_ir_x), ]
pcair <- cbind(pca_ir_x, meta_r_ord[match(rownames(pca_ir_x), meta_r_ord$`Sample Name_RNA`), ])

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = ID, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors_ir) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = HSCT_responder, 
                label = paste(Time, ID, sep = "_"))) +
  ggtitle("PCA RNAseq biopsies") + 
  guides(col = guide_legend(title = "Responders")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%"))

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = CD_Aftected_area, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Afected area")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")


ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = Exact_location, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Region")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = Involved_Healthy, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Involved area")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")


# PCA TRIM (without barcelona)
barcelona <- grepl("w", colnames(expr))
counts_woBarcelona <- expr[, !barcelona]
counts_woBarcelona <- counts_woBarcelona[rowSums(counts_woBarcelona) != 0, ]
pca_ir <- prcomp(t(counts_woBarcelona), scale. = TRUE)
pca_ir_x <- as.data.frame(pca_ir$x)
pca_ir_var <- round(summary(pca_ir)$importance[2, ]*100, digits = 2)

pca_ir_x <- pca_ir_x[rownames(pca_ir_x) %in% meta_r$`Sample Name_RNA`, ]
meta_r_ord <- meta_r[meta_r$`Sample Name_RNA` %in% rownames(pca_ir_x), ]
pcair <- cbind(pca_ir_x, meta_r_ord[match(rownames(pca_ir_x), meta_r_ord$`Sample Name_RNA`), ])

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = ID, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors_ir) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = HSCT_responder, 
                label = paste(Time, ID, sep = "_"))) +
  ggtitle("PCA RNAseq biopsies") + 
  guides(col = guide_legend(title = "Responders")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%"))

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = CD_Aftected_area, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")

ggplot(pcair) +
  geom_text(aes(PC1, PC2, col = Involved_Healthy, label = paste(ID, Time, sep = "_"))) + 
  guides(col = guide_legend(title="Patient")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(paste("PC1", pca_ir_var[1], "%")) +
  ylab(paste("PC2", pca_ir_var[2], "%")) + 
  ggtitle("PCA RNAseq biopsies")

dev.off()
