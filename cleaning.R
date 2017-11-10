intestinal <- "intestinal_16S"
stool <- "stools_16S"
rna <- "intestinal_RNAseq"
source("helper_functions.R")

# Read the intestinal otus table
otus_table_i <- read.csv(file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
                         stringsAsFactors = FALSE, row.names = 1, 
                         check.names = FALSE)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

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
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(file_meta_i, row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE)
file_meta_r <- file.path(rna, "metadata.csv")
meta_r <- read.csv(file_meta_r, row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE)

# Reorder by patient and time
meta_s <- meta_s[order(meta_s$Patient_ID, meta_s$Time), ]
meta_i <- meta_i[order(meta_i$Patient_ID, meta_i$Time), ]
meta_r <- meta_r[order(meta_r$Patient_ID, meta_r$Time), ]

# Find common patients
comPatient <- intersect(meta_i$Patient_ID, meta_s$Patient_ID)
# Controls with same name but they are different people
comPatient <- comPatient[!grepl("^C", comPatient)] 
comTime <- intersect(meta_i$Time, meta_s$Time)

# Keep only the common patients and times
com_meta_i <- meta_i[meta_i$Patient_ID %in% comPatient & meta_i$Time %in% comTime, ]
com_meta_s <- meta_s[meta_s$Patient_ID %in% comPatient & meta_s$Time %in% comTime, ]

# Delete rows of patients which are not in common between the datasets
keep_i <- rownames(meta_i)[meta_i$Patient_ID %in% comPatient & 
                             meta_i$Time %in% comTime]
keep_s <- rownames(meta_s)[meta_s$Patient_ID %in% comPatient & 
                             meta_s$Time %in% comTime]
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
add_s <- meltQ[meltQ$value > 0 & meltL$value,]
remove_i <- meltQ[meltQ$value > 0 & !meltL$value,]

# Remove the stools that should be removed because we don't have intestinal 
# samples
rownames(remove_s) <- seq_len(nrow(remove_s))

remove_samples_s <- apply(remove_s, 1, function(x){
  rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]]
})

com_otus_table_s <- com_otus_table_s[, 
                                     !colnames(com_otus_table_s) %in% remove_samples_s]

# Add stool samples as much as need
rownames(add_s) <- seq_len(nrow(add_s))

add_samples_s <- apply(add_s, 1, function(x){
  rep(rownames(meta_s)[meta_s$Patient_ID %in% x[1] & meta_s$Time %in% x[2]], x[3])
})

com_otus_table_s <- com_otus_table_s[, 
                                     c(colnames(com_otus_table_s), unlist(add_samples_s))]

# Remove intestinal data
rownames(remove_i) <- seq_len(nrow(remove_i))

remove_samples_i <- apply(remove_i, 1, function(x){
  rownames(meta_i)[meta_i$Patient_ID %in% x[1] & meta_i$Time %in% x[2]]
})

com_otus_table_i <- com_otus_table_i[, 
                                 !colnames(com_otus_table_i) %in% unlist(remove_samples_i)]

# Reorder so that the stools samples and the intestinal samples are in the 
# same order including Time and Patient_ID

char_com_i <- meta_i[colnames(com_otus_table_i), c("Patient_ID", "Time")]

names_clean <- gsub("(.+)\\.[0-9]$", "\\1", colnames(com_otus_table_s))
char_com_s <- meta_s[names_clean, c("Patient_ID", "Time")]
o <- order(char_com_s$Patient_ID, char_com_i$Time)

otus_s <- com_otus_table_s[, o]
otus_i <- com_otus_table_i

# We transpose the data because it requires the data in column for variable, row
# for sample and we remove thosw which are all empty
keep_otus_s <- apply(t(otus_s), 2, sd) != 0
keep_otus_i <- apply(t(otus_i), 2, sd) != 0
otus_s_f <- t(otus_s)[, keep_otus_s]
otus_i_f <- t(otus_i)[, keep_otus_i]

# Clean the metadata
meta <- com_meta_i[colnames(otus_i),]
meta$HSCT_responder[meta$HSCT_responder == "C"] <- NA
meta$Active_area[meta$Active_area == ""] <- NA

# Remove non informative variables
meta <- meta[, apply(meta, 2, function(x){length(unique(x)) != 1})]


## Find the otus that are equivalent between datasets
comb <- expand.grid(rownames(otus_tax_i[keep_otus_i, ]), 
                    rownames(otus_tax_s[keep_otus_s, ]), stringsAsFactors = FALSE)
colnames(comb) <- c("intestinal", "stools")

eq <- apply(comb, 1, function(z){
  y <- z[2]
  x <- z[1]
  # If there is any NA then they are nor precise enough to say they are the same
  # OTU
  all(otus_tax_i[keep_otus_i, ][x, ] == otus_tax_s[keep_otus_s, ][y, ]) 
})

eqOTUS <- comb[eq & !is.na(eq), ]
eqOTUS <- droplevels(eqOTUS)
rownames(eqOTUS) <- seq_len(nrow(eqOTUS))

# Write the files
write.csv(otus_s_f, row.names = FALSE, 
          file = "stools_16S/otus_coherent.csv")
write.csv(otus_i_f, row.names = FALSE, 
          file = "intestinal_16S//otus_coherent.csv")
write.csv(meta, file = "meta_coherent.csv")
write.csv(otus_tax_i[keep_otus_i, ], 
          file = file.path(intestinal, "taxonomy.csv"), row.names = TRUE)
write.csv(otus_tax_s[keep_otus_s, ], 
          file = file.path(stool, "taxonomy.csv"), row.names = TRUE)
write.csv(eqOTUS, "equivalent_otus.csv", row.names = FALSE)