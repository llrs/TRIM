library("integration")
library("SummarizedExperiment")
library("reshape")

intestinal <- "intestinal_16S"
stool <- "stools_16S"
rna <- "intestinal_RNAseq"
today <- format(Sys.time(), "%Y%m%d")

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Read the stools OTUs
otus_table_s <- read.delim(
  file.path(stool, "OTUs-Table-refined-stools.tab"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_s <- otus_table_s[, ncol(otus_table_s)]
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

# Extract the taxonomy and format it properly
otus_tax_s <- taxonomy(tax_s, rownames(otus_table_s))

# Load the input data
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(
  file_meta_s, check.names = FALSE, row.names = 1,
  stringsAsFactors = FALSE
)
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_complete.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

# Correct the swapped samples and match metadata
expr <- norm_expr_colnames(expr)
pdf(paste0("Figures/", today, "_quality.pdf"))
counts <- colSums(otus_table_i)
counts_ord <- counts[order(counts)]
code <- meta_i$Sample_Code
names(code) <- rownames(meta_i)
barplot(counts_ord, col = ifelse(grepl("w", code[names(counts_ord)]), "red", "black"),
        main = "Total counts otus biopsies", xlab = "Samples", ylab = "counts",
        border = NA, 
        sub = paste(round(sum(counts_ord < 10000)/length(counts_ord)*100), "%"))
abline(h = 10000, col = "green")
abline(h = mean(counts), col = "orange")
abline(h = 25000)

o <- apply(otus_table_i, 2, function(x) sum(x != 0))
mean(o)
p <- data.frame(OTUs = o, Abundance = counts)
ggplot(p, aes(Abundance, OTUs)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Biopsies") +
  xlim(c(0, max(p$Abundance)))

p2 <- cbind(p, meta_i[match(rownames(p), rownames(meta_i)), ])
# Remove a genus because it is an outlier
ggplot(p2[-88, ], aes(Abundance, OTUs)) +
  geom_point(aes(col = Time)) +
  geom_smooth() +
  ggtitle("Biopsies", sub = "Removing outlier sample") +
  xlim(c(0, max(p$Abundance[-88])))

counts <- colSums(otus_table_s)
counts_ord <- counts[order(counts)]
barplot(counts_ord, col = "black",
        main = "Total counts otus stools", xlab = "Samples", ylab = "counts",
        border = NA)
abline(h = 100000, col = "green")
abline(h = mean(counts), col = "orange")
abline(h = 50000)

o <- apply(otus_table_s, 2, function(x) sum(x != 0))
p <- data.frame(OTUs = o, Abundance = counts)
ggplot(p, aes(Abundance, OTUs)) +
  geom_point() +
  geom_smooth() +
  ggtitle("Stools") +
  xlim(c(0, max(p$Abundance)))

counts <- colSums(expr)
counts_ord <- counts[order(counts)]
barplot(counts_ord, col = ifelse(grepl("w", names(counts)), "red", "black"),
        main = "Total counts biopsies RNAseq", xlab = "Samples",
        border = NA, ylab = "counts")
abline(h = 100000, col = "green")
abline(h = mean(counts), col = "orange")
abline(h = 50000)

o <- apply(expr, 2, function(x) sum(x != 0))
p <- data.frame(Genes = o, Abundance = counts)
ggplot(p, aes(Abundance, Genes)) +
  geom_point() +
  geom_smooth() +
  ggtitle("RNA")

dev.off()

# Reorder by patient and time
meta_s <- meta_s[order(meta_s$Patient_ID, meta_s$Time), ]
meta_i <- meta_i[order(meta_i$Patient_ID, meta_i$Time), ]
meta_r <- meta_r[order(meta_r$Patient_ID, meta_r$Time), ]

# Correct samples metadata
meta_r <- meta_r_norm(meta_r)
meta_i <- meta_i_norm(meta_i)

# Find common patients
comPatient <- intersect(meta_i$Patient_ID, meta_s$Patient_ID)
# Controls with same name but they are different people
comPatient <- comPatient[!grepl("^C", comPatient)]
comTime <- intersect(meta_i$Time, meta_s$Time)

# Keep only the common patients and times
com_meta_i <- meta_i[meta_i$Patient_ID %in% comPatient & meta_i$Time %in% comTime, ]
com_meta_s <- meta_s[meta_s$Patient_ID %in% comPatient & meta_s$Time %in% comTime, ]

comPatient <- intersect(com_meta_i$Patient_ID, com_meta_s$Patient_ID)
com_meta_i <- com_meta_i[com_meta_i$Patient_ID %in% comPatient, ]
com_meta_s <- com_meta_s[com_meta_s$Patient_ID %in% comPatient, ]

# Delete rows of patients which are not in common between the datasets
keep_i <- rownames(meta_i)[meta_i$Patient_ID %in% comPatient &
  meta_i$Time %in% comTime]
keep_s <- rownames(meta_s)[meta_s$Patient_ID %in% comPatient &
  meta_s$Time %in% comTime]

keep_i <- keep_i[keep_i %in% colnames(otus_table_i)]
keep_s <- keep_s[keep_s %in% colnames(otus_table_s)]

com_otus_table_i <- otus_table_i[, keep_i]
com_otus_table_s <- otus_table_s[, keep_s]

# Samples per patient and Time
tab_s <- table(com_meta_s$Patient_ID, com_meta_s$Time)
tab_i <- table(com_meta_i$Patient_ID, com_meta_i$Time)

# Compare if we have more intestinal samples than stools for the same time and
# patient
moreS <- tab_i > tab_s & tab_s != 0
QmoreS <- tab_i - tab_s

# Extract the ids of the data we need to add or remove
meltQ <- melt(QmoreS, varnames = c("Patient_ID", "Time"))
meltL <- melt(moreS, varnames = c("Patient_ID", "Time"))
remove_s <- meltQ[meltQ$value < 0, ]
add_s <- meltQ[meltQ$value > 0 & meltL$value, ]
remove_i <- meltQ[meltQ$value > 0 & !meltL$value, ]

# Remove the stools that should be removed because we don't have intestinal
# samples
rownames(remove_s) <- seq_len(nrow(remove_s))

remove_samples_s <- apply(remove_s, 1, function(x) {
  rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]]
})

com_otus_table_s <- com_otus_table_s[
  ,
  !colnames(com_otus_table_s) %in% remove_samples_s
]

# Add stool samples as much as need
rownames(add_s) <- seq_len(nrow(add_s))

add_samples_s <- apply(add_s, 1, function(x) {
  rep(rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]], 
      x[3])
})

com_otus_table_s <- com_otus_table_s[
  , c(colnames(com_otus_table_s), unlist(add_samples_s))]

# Remove intestinal data
rownames(remove_i) <- seq_len(nrow(remove_i))

remove_samples_i <- apply(remove_i, 1, function(x) {
  rownames(meta_i)[meta_i$Patient_ID %in% x[1] & meta_i$Time %in% x[2]]
})

com_otus_table_i <- com_otus_table_i[
  , !colnames(com_otus_table_i) %in% unlist(remove_samples_i)]

# Reorder so that the stools samples and the intestinal samples are in the
# same order including Time and Patient_ID

char_com_i <- meta_i[colnames(com_otus_table_i), c("Patient_ID", "Time")]

names_clean <- gsub("(.+)\\.[0-9]$", "\\1", colnames(com_otus_table_s))
char_com_s <- meta_s[names_clean, c("Patient_ID", "Time")]
o <- order(char_com_s$Patient_ID, char_com_i$Time)

otus_s <- com_otus_table_s[, o]
otus_i <- com_otus_table_i

# We transpose the data because it requires the data in column for variable, row
# for sample and we remove those which are all empty
keep_otus_s <- apply(t(otus_s), 2, sd) != 0
keep_otus_i <- apply(t(otus_i), 2, sd) != 0
otus_s_f <- t(otus_s)[, keep_otus_s]
otus_i_f <- t(otus_i)[, keep_otus_i]
tax_i <- otus_tax_i[keep_otus_i, ]
tax_s <- otus_tax_s[keep_otus_s, ]

# Clean the metadata
meta <- com_meta_i[colnames(otus_i), ]
meta$HSCT_responder[meta$HSCT_responder == "C"] <- NA
meta$Active_area[meta$Active_area == ""] <- NA

# Remove non informative variables
meta <- meta[, apply(meta, 2, function(x) {
  length(unique(x)) != 1
})]
meta$Active_area[meta$Active_area == ""] <- NA
meta$ID <- meta$Patient_ID
meta$ID[meta$Patient_ID %in% c("15", "23")] <- "15/23"
meta$ID[meta$Patient_ID %in% c("33", "36")] <- "33/36"
meta$ID[meta$Patient_ID %in% c("29", "35")] <- "29/35"
meta$ID <- as.factor(meta$ID)

# Pre transplant
meta$Transplant <- "Post" #
meta$Transplant[meta$Patient_ID %in% c("15", "33", "29")] <- "Pre"

## Find the otus that are equivalent between datasets
comb <- expand.grid(
  rownames(otus_tax_i[keep_otus_i, ]),
  rownames(otus_tax_s[keep_otus_s, ]), stringsAsFactors = FALSE
)
colnames(comb) <- c("intestinal", "stools")

eq <- apply(comb, 1, function(z) {
  y <- z[2]
  x <- z[1]
  # If there is any NA then they are nor precise enough to say they are the same
  # OTU
  sum(otus_tax_i[keep_otus_i, ][x, ] == otus_tax_s[keep_otus_s, ][y, ])
})

eqOTUS <- comb[eq >= 7 & !is.na(eq), ]
eqOTUS <- droplevels(eqOTUS)
rownames(eqOTUS) <- seq_len(nrow(eqOTUS))

eqGenera <- comb[eq >= 6 & !is.na(eq), ]
eqGenera <- droplevels(eqGenera)
rownames(eqGenera) <- seq_len(nrow(eqGenera))

# Create the SummarizedExperiment objects
meta_i <- meta_i[match(colnames(otus_table_i), rownames(meta_i)), ]
SE_i <- SummarizedExperiment(
  assays = SimpleList(otus = as.matrix(otus_table_i)),
  colData = meta_i,
  rowData = otus_tax_i
)

meta_s <- meta_s[match(colnames(otus_table_s), rownames(meta_s)), ]
SE_s <- SummarizedExperiment(
  assays = SimpleList(otus = as.matrix(otus_table_s)),
  colData = meta_s,
  rowData = otus_tax_s
)
meta_r <- meta_r[match(colnames(expr), meta_r$`Sample Name_RNA`), ]
rownames(meta_r) <- NULL
SE_expr <- SummarizedExperiment(
  assays = SimpleList(expr = as.matrix(expr)),
  colData = meta_r
)

# prepMultiAssay(list(SE_i, SE_s, SE_expr), meta_r, )

list_SE <- list(
  "intestinal_16S" = SE_i,
  "stools_16S" = SE_s,
  "intestinal_RNAseq" = SE_expr
)
assay <- rep(names(list_SE), sapply(list_SE, function(x) {
  ncol(assay(x))
}))
colnames <- unlist(sapply(list_SE, function(x) {
  colnames(assay(x))
}), use.names = FALSE)
sampleMap <- data.frame(
  assay = as.factor(assay), primary = NA,
  colname = colnames
)
# sampleMap$primary[sampleMap$assay == sampleMap$assay[1]]
# <- meta_r[, meta_r$`Sample Name_RNA` == sampleMap$colname]
#

# Find the samples that we have microbiota and expression
i_names <- rownames(meta_i)[meta_i$Sample_Code %in% meta_r$Sample_Code_uDNA]
i_IDs <- meta_i$Sample_Code[rownames(meta_i) %in% i_names]
e_IDs <- meta_r$`Sample Name_RNA`[meta_r$Sample_Code_uDNA %in% i_IDs]
o_i <- otus_table_i[, i_names]
e_i <- expr[, e_IDs]

# Write the files
write.csv(
  otus_s_f, row.names = FALSE,
  file = "stools_16S/otus_coherent.csv"
)
write.csv(
  otus_i_f, row.names = FALSE,
  file = "intestinal_16S//otus_coherent.csv"
)
write.csv(meta, file = "meta_coherent.csv")
write.csv(tax_i, file = file.path(intestinal, "taxonomy.csv"), row.names = TRUE)
write.csv(tax_s, file = file.path(stool, "taxonomy.csv"), row.names = TRUE)
write.csv(eqOTUS, "equivalent_otus.csv", row.names = FALSE)
write.csv(eqGenera, "equivalent_genus.csv", row.names = FALSE)
