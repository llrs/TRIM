# Summarize by taxa
library("metagenomeSeq")
library("integration")
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("forcats")

wd <- setwd("..")
today <- format(Sys.time(), "%Y%m%d")

intestinal <- "intestinal_16S"
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

# Read the metadata for each type of sample
rna <- "intestinal_RNAseq"
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(wd)

# Clean the metadata
meta_r <- meta_r_norm(meta_r)

# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Subset and reorder
meta_r <- meta_r[match(colnames(otus_table_i), meta_r$Seq_code_uDNA), ]
rownames(meta_r) <- meta_r$Seq_code_uDNA

# Convert to genus 
MR_i <- newMRexperiment(
  otus_table_i,
  phenoData = AnnotatedDataFrame(meta_r),
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i))
)

genus <- aggTax(MR_i, log = FALSE, lvl = "Genus", norm = FALSE, out = "matrix")


out <- genus %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column() %>% 
  gather(... = -"rowname") %>% 
  rename(Genus = rowname,
         Sample = key,
         Count = value) %>% 
  group_by(Sample) %>% 
  filter(Genus != "no_match") %>% 
  mutate(Abundance = Count/sum(Count)*100) %>% 
  ungroup()

topMicro <- out %>% 
  group_by(Genus) %>% 
  summarise(mA=mean(Abundance)) %>% 
  arrange(desc(mA)) %>% 
  top_n(10, mA)

m <- out[out$Genus %in% topMicro$Genus, ]
o <- merge(m, meta_r, by.x = "Sample", by.y = "Seq_code_uDNA", all.x = TRUE)

o %>% 
  filter(IBD != "CONTROL") %>% 
  mutate(Time = case_when(Time == "S0" ~ "T0", TRUE ~ Time)) %>% 
  droplevels() %>% 
  ggplot() +
  geom_tile(aes(lvls_reorder(Time, c(1, 3, 4, 5, 6, 7, 8, 2)), Genus, fill = Abundance)) +
  facet_wrap(~ID, nrow = 2) +
  scale_fill_gradient2() +
  labs(fill = "Abundance (%)", x = "Time") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  