intestinal <- "intestinal_16S"
stool <- "stools_16S"
rna <- "intestinal_RNAseq"

today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("grid")

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Read the stools OTUs
otus_table_s <- read.delim(
  file.path(stool, "OTUs-Table-refined-stools.tab"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
otus_table_s <- otus_table_s[, -ncol(otus_table_s)]

# Load the RNAseq
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

# Correct the swapped samples
position <- c(grep("33-T52-TTR-CIA", colnames(expr)), 
              grep("33-T52-TTR-IIA", colnames(expr)))
colnames(expr)[position] <- colnames(expr)[rev(position)]

colnames(expr) <- toupper(colnames(expr))

# Read the metadata for each type of sample
file_meta_s <- "stools_16S/db_stool_samples_microbiome_abstract_RUN3def.txt"
meta_s <- read.delim(
  file_meta_s, check.names = FALSE, row.names = 1,
  stringsAsFactors = FALSE
)
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)

file_meta_r <- file.path(rna, "metadata_28032018.csv")
meta_r <- read.table(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, sep = ";",
  na.strings = c("NA", ""),
)
colnames(meta_r) <- meta_r[1, ]
meta_r <- meta_r[-1, ]

# Clean the metadata
meta_i <- meta_i_norm(meta_i)
meta_s <- meta_s_norm(meta_s)
meta_r <- meta_r_norm(meta_r)

# Check that the names of the matrices and the names of the metadata are 
# equivalent
if (nrow(meta_i[colnames(otus_table_i), ]) != nrow(meta_i)){
  meta_i <- meta_i[colnames(otus_table_i), ]
}

if (nrow(meta_s[colnames(otus_table_s), ]) != nrow(meta_s)){
  meta_s <- meta_s[colnames(otus_table_s), ]
}

if (nrow(meta_r[colnames(expr), ]) != nrow(meta_r)){
  meta_r <- meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ]
  if (any(!colnames(expr) %in% meta_r$`Sample Name_RNA`)){
    missing <- colnames(expr)[!colnames(expr) %in% meta_r$`Sample Name_RNA`]
    if (all(grepl("-w[0-9]+", missing, ignore.case = TRUE))){
      missing <- missing[!grepl("-w[0-9]+", missing, ignore.case = TRUE)]
      message("Barcelona samples")
    } else {
      stop("Missing metadata")
    }
  }
}

# First by patient
i_pat <- unique(as.character(meta_i$ID))
s_pat <- unique(as.character(meta_s$ID))
r_pat <- unique(as.character(meta_r$ID))
patient_name <- unique(c(i_pat, s_pat, r_pat))
zeroes <- rep(0, length(patient_name))
df <- sapply(patient_name, function(x){
  Biopsies_RNAseq <- as.numeric(x %in% meta_r$ID)
  Biopsies_16S <- as.numeric(x %in% meta_i$ID)
  Stools_16S <- as.numeric(x %in% meta_s$ID)
  c("Biopsies_RNAseq" = Biopsies_RNAseq, "Biopsies_16S" = Biopsies_16S, 
    "Stools_16S" = Stools_16S)
})

df <- as.data.frame(t(df))

# Controls from biopsies are not the same as controls from biopsies
df <- cbind(df, Controls = as.numeric(grepl("^C|[A-Z]", rownames(df))))
df$Controls[grepl("^C|[A-Z]", rownames(df))] <- 1
subdf <- df[df$Controls == 1, ]
df[df$Controls == 1, "Stools_16S"] <- 0 
stool_sub <- subdf$Stools_16S
subdf[,] <- 0
subdf$Stools_16S <- stool_sub
subdf <- subdf[rowSums(subdf) != 0, ]
subdf$Controls <- 1
rownames(subdf) <- paste0(rownames(subdf), "_2")
df <- rbind(df, subdf)
df <- df[!grepl("^(JR|NP)$", rownames(df)), ]

metadata <- data.frame(sets = colnames(df), "Type" = c("Biopsies", "Biopsies", "Stools", "Depends"))
library("UpSetR")
pdf("Figures/patients_dataset.pdf")
upset(df, order.by = "freq", line.size = NA, sets.x.label = "Patients samples",
      mainbar.y.label = "Patients intersection", 
      set.metadata = list(data = metadata, plots = list(list(type = "matrix_rows", 
                                                             column = "Type", 
                                                             colors = c(Biopsies = "green", Stools = "blue", Depends = "Grey"), 
                                                alpha = 0.25))))
dev.off()

# Then by biopsies
# meta_i$Sample_code ==  meta_r$Sample_Code_uDNA

# As the name of nam is in the same order than meta_r we can just iterate 
# through meta_r
incorporate <- function(vector){
  vector <- as.character(vector)
  v <- unique(vector)
  v <- v[!is.na(v)]
  m2 <- sapply(v, function(x){as.numeric(vector %in% x)})
  as.data.frame(m2)
}

m_time <- incorporate(meta_r$Time)
time <- colnames(m_time)

m_location <- incorporate(meta_r$Exact_location)
location <- colnames(m_location)
m <- cbind(m_time, m_location)

st <- data.frame("Biopsies_16S" = as.numeric(!is.na(meta_r$Sample_Code_uDNA)), 
                 "Biopsies_RNAseq" = as.numeric(!is.na(meta_r$`Sample Name_RNA`)))
m <- cbind(m, st)
st <- colnames(st)

m2 <- cbind(m, meta_r[, c("Time", "Involved_Healthy", "SEX", "Patient_ID", "Sample_Code_uDNA")])
m_controls <- m2[m2$C == 1, ]
m_IBD <- m2[m2$C == 0, ]

# Reorder 
time <- time[-match("C", time)]
time <- time[c(1, 8, 3, 7, 2, 5, 6)]

location <- location[c(1, 4, 6, 5, 2, 3)]


pdf("Figures/samples_dataset.pdf")
upset(m_IBD, order.by = "freq", sets = st, line.size = NA)
grid.text("TRIM", x = 0.65, y=0.95, gp=gpar(fontsize=20))
upset(m_IBD, order.by = "freq", sets = location, line.size = NA, keep.order = TRUE)
grid.text("Localization", x = 0.65, y=0.95, gp=gpar(fontsize=20))
upset(m_IBD, order.by = "freq", sets = time, line.size = NA, keep.order = TRUE)
grid.text("Time", x = 0.65, y=0.95, gp=gpar(fontsize=20))
upset(m_IBD, order.by = "freq", sets = c(st, location), line.size = NA)
upset(m_IBD, order.by = "freq", sets = c(st, time), line.size = NA)
upset(m_IBD, order.by = "freq", sets = c(st, location, time), line.size = NA)

upset(m_controls, sets = st, line.size = NA)
grid.text("Samples", x = 0.65, y = 0.95, gp=gpar(fontsize=15))
keep <- sapply(location, function(x){any(m_controls[, x] != 0)})
upset(m_controls, order.by = "freq", sets = location[keep], line.size = NA, keep.order = TRUE)
grid.text("Localization", x = 0.65, y=0.95, gp=gpar(fontsize=20))
grid.text("(Controls)", x = 0.79, y = 0.95, gp=gpar(fontsize=15))
upset(m_controls, order.by = "freq", sets = c(st, location[keep]), line.size = NA)
grid.text("Controls", x = 0.670, y = 0.95, gp=gpar(fontsize=15))
dev.off()

# Calculate prevalence plots
microab <- read.csv("AllSamples_Time_MicrobesAbund.csv", stringsAsFactors = FALSE)
mm_type <- incorporate(microab$Sample_Type)
mm_time <- incorporate(microab$Time)
mm_patient <- incorporate(microab$Patient_ID)
mm_response <- incorporate(microab$HSCT_responder)
mm <- cbind(mm_type, mm_response, mm_time, 
      "g__Megasphaera" = ifelse(microab$g__Megasphaera > 0.5, 1, 0),
      "g__Streptococcus" = ifelse(microab$g__Streptococcus > 0.5, 1, 0),
      "s__prausnitzii" = ifelse(microab$s__prausnitzii > 0.5, 1, 0))
colnames(mm)[match("B", colnames(mm))] <- "Biopsy"
colnames(mm)[match("S", colnames(mm))] <- "Stools"
colnames(mm)[match("YES", colnames(mm))] <- "Responder"
mm <- mm[, -match(c("NO"), colnames(mm))]
sets <- colnames(mm)
mm <- cbind(mm, mm_patient, microab)

# Extract from where is each sample
subs <- strsplit(mm$Sample.code..TOTAL., "_")
tab <- sapply(subs, function(x){
  if (length(x) <5){
    0
  } else {
    substr(x[5], 1, 1)
  }
})
mm_area <- data.frame("Colon" = ifelse(tab == "C", 1, 0),
                      "Ileum" = ifelse(tab == "I", 1, 0))
mm <- cbind(mm, mm_area)
library("data.table")
setDT(mm)
cols <- c("Stools", "Responder", "T0", "T26", "T52",
          "g__Megasphaera", "g__Streptococcus", "s__prausnitzii", "Colon", "Ileum")
count <- function(x){
  ifelse(x != 0, 1, 0)
}
# Summarize by column?
# mm[, .SD[, .N], by = c("Patient_ID", "Sample_Type", "Stools", "Colon", "Ileum"), .SDcols = cols]

d <- list()
for (ID in as.character(unique(mm$Patient_ID))){
   x <- colSums(mm[mm$Patient_ID == ID & mm$Sample_Type == "B", 
                   c("Colon", "Ileum", "g__Megasphaera", "g__Streptococcus", "s__prausnitzii")], 
                na.rm = TRUE)
   d[[ID]] <- ifelse(x != 0, 1, 0)
}
m_samples_biopsies <- t(simplify2array(d))
m_samples_biopsies <- m_samples_biopsies[rowSums(m_samples_biopsies) != 0, ]
m_samples_biopsies <- as.data.frame(m_samples_biopsies)


pdf("Figures/abundance.pdf")
# upset(mm[mm$Biopsy == 1, ], 
#       sets = c("Responder", "g__Megasphaera", "g__Streptococcus", 
#                "s__prausnitzii"), order.by = "freq", 
#       line.size = NA)
# grid.text("Biopsies", x = 0.65, y=0.95, gp=gpar(fontsize=20))
# upset(mm[mm$Biopsy == 0, ], 
#       sets = c("Responder", "g__Megasphaera", "g__Streptococcus", 
#                "s__prausnitzii"), order.by = "freq", 
#       line.size = NA)
# grid.text("Stools", x = 0.65, y=0.95, gp=gpar(fontsize=20))
# upset(mm, 
#       sets = c("Biopsy", "Responder", "g__Megasphaera", "g__Streptococcus", 
#                "s__prausnitzii"), order.by = "freq", 
#       line.size = NA)
upset(m_samples_biopsies, order.by = "freq", 
      sets = c("Colon", "Ileum", "s__prausnitzii"), 
      line.size = NA)
grid.text("Prausnitzii", x = 0.65, y=0.95, gp=gpar(fontsize=20))
upset(m_samples_biopsies, 
      sets = c("Colon", "Ileum", "g__Megasphaera"), order.by = "freq", 
      line.size = NA)
grid.text("Megasphaera", x = 0.65, y=0.95, gp=gpar(fontsize=20))

upset(m_samples_biopsies, 
      sets = c("Colon", "Ileum", "g__Streptococcus"), order.by = "freq", 
      line.size = NA)
grid.text("Streptococcus", x = 0.65, y=0.95, gp=gpar(fontsize=20))


d <- list()
for (ID in as.character(unique(mm$Patient_ID))){
  x <- colSums(mm[mm$Patient_ID == ID & mm$Sample_Type == "S", 
                  c("g__Megasphaera", "g__Streptococcus", "s__prausnitzii")], 
               na.rm = TRUE)
  d[[ID]] <- ifelse(x != 0, 1, 0)
}
m_samples_stools <- t(simplify2array(d))
m_samples_stools <- m_samples_stools[rowSums(m_samples_stools) != 0, ]
m_samples_stools <- as.data.frame(m_samples_stools)

upset(m_samples_stools, 
      sets = colnames(m_samples_stools), order.by = "freq", 
      line.size = NA)
grid.text("Taxas in stools", x = 0.65, y=0.95, gp=gpar(fontsize=20))



# upset(mm[, c("Stools", "Colon", "Ileum", "g__Megasphaera", "g__Streptococcus", 
#              "s__prausnitzii")], 
#       nintersects = NA, order.by = "freq",
#       line.size = NA)
# upset(mm, sets = sets, order.by = "freq", line.size = NA)
# grid.text("Time influence", x = 0.65, y=0.95, gp=gpar(fontsize=20))
# upset(mm, sets = c(unique(microab$Patient_ID), sets), order.by = c("degree"), 
#       nintersects = NA, line.size = NA,
#       queries = list(list(query = elements, params = list("Patient_ID"), 
#                         color = "blue", active = TRUE)))
# grid.text("By patient", x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()