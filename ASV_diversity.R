library("integration")
library("phyloseq")
library("dplyr")
library("tidyr")
library("ggplot2")

rna <- "intestinal_RNAseq"

file_meta <- file.path(rna, "metadata_25042018.csv")
meta <- read.delim(
  file_meta, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)
meta <- meta_r_norm(meta)
meta_c <- meta[meta$Time == "C", ]

seqtab.nochim <- readRDS("intestinal_16S_RNAseq_metadb/dada2.RDS")
ASV_counts <- t(seqtab.nochim)
rownames(ASV_counts) <- NULL

df <- strcapture("filt_([0-9]+)\\.(.+)@.\\.fastq$",
           colnames(ASV_counts), 
           proto = data.frame(code = numeric(), sample = character()))
df$file <- colnames(ASV_counts)


df_c <- df[df$sample %in% meta_c$Seq_code_uDNA, ]

meta_c$Sample_Code_uDNA <- gsub("_", "-", meta_c$Sample_Code_uDNA)

samples_aTNF <- c("C1-T-DM-SIG", "C10-T-DM-ILI", "C10-T-DM-SIG", "C3-T-DM-SIG", 
                  "C4-T-DM-SIG", "C6-T-DM-ILI", "C6-T-DM-SIG", "C7-T-DM-ILI", "C7-T-DM-SIG", 
                  "C8-T-DM-ILI", "C8-T-DM-SIG", "C9-T-DM-ILI", "C9-T-DM-SIG")
m <- meta_c[meta_c$Sample_Code_uDNA %in% samples_aTNF, ]
ctrls_aTNF <- df_c[df_c$sample %in% m$Seq_code_uDNA, ]
df_meta <- merge(m, df_c[, 1:2], by.x = "Seq_code_uDNA", by.y = "sample")

saveRDS(df_meta, "meta_shared_controls.RDS")

asv_samples <- ASV_counts[, ctrls_aTNF$file]
colnames(asv_samples) <- ctrls_aTNF$sample[match(ctrls_aTNF$file, colnames(asv_samples))]
m <- m[match(colnames(asv_samples), m$Seq_code_uDNA), ]
stopifnot(m$Seq_code_uDNA == colnames(asv_samples))
rownames(m) <- colnames(asv_samples)
phyloseq <- phyloseq(otu_table(asv_samples, taxa_are_rows = TRUE),
                     sample_data(m)
                     # tax_table(as.matrix(family))
)
theme_set(theme_bw())
alpha_meas <- c("Simpson", "Shannon")
richness <- estimate_richness(phyloseq)
richness <- cbind(richness, m)
r2 <- pivot_longer(richness, cols = Chao1:Fisher, names_to = "Alpha diversity")
richness_rel <- filter(r2, `Alpha diversity` %in% c("Shannon", "Simpson")) %>%
  mutate(effective = case_when(
    `Alpha diversity` == "Shannon" ~ exp(value),
    `Alpha diversity` == "Simpson" ~ 1/value,
  ))
saveRDS(richness_rel, "alpha_diversity_aTNF_controls.RDS")
saveRDS(richness_rel, "_aTNF_controls.RDS")


# Diversity of all samples
ASV_counts <- ASV_counts[, colnames(ASV_counts) %in% df$file]
colnames(ASV_counts) <- df$sample[match(df$file, colnames(ASV_counts))]
meta_all <- meta[match(colnames(ASV_counts), meta$Seq_code_uDNA), ]
rownames(meta_all) <- colnames(ASV_counts)
phyloseq <- phyloseq(otu_table(ASV_counts, taxa_are_rows = TRUE),
                     sample_data(meta_all)
                     # tax_table(as.matrix(family))
)
richness <- estimate_richness(phyloseq)
richness <- cbind(richness, meta_all)
r2 <- pivot_longer(richness, cols = Chao1:Fisher, names_to = "Alpha diversity")
richness_rel <- filter(r2, `Alpha diversity` %in% c("Shannon", "Simpson")) %>%
  mutate(effective = case_when(
    `Alpha diversity` == "Shannon" ~ exp(value),
    `Alpha diversity` == "Simpson" ~ 1/value,
  ),
  ileum = ifelse(Exact_location %in% "ILEUM", "ileum", "colon"),
  IBD = forcats::fct_relevel(IBD, "CONTROL", "CD", "UC"))
ggplot(richness_rel, aes(IBD, effective)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_jitter(height = 0) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank()) +
  theme_minimal()
ggsave("Figures/alpha_diversity.png")
# To filter just for our samples
meta_r <- readRDS("intestinal_16S_RNAseq_metadb/meta.RDS")
richness_rel %>% 
  filter(Seq_code_uDNA %in% meta_r$Seq_code_uDNA) %>% dim()
  ggplot(aes(IBD, effective)) +
  geom_boxplot(alpha = 0, outlier.size = 0) +
  geom_jitter(height = 0) +
  facet_grid(`Alpha diversity` ~ ileum, scales = "free", drop = TRUE) +
  labs(y = "Alpha diversity", x = element_blank()) +
  theme_minimal()
ggsave("Figures/alpha_diversity.png")


aTNF_ASV <- readRDS("../design_ngs/data_out/ASV_sequences.RDS")
TRIM_ASV <- colnames(seqtab.nochim) 
hist(nchar(TRIM_ASV), breaks = 20, xlim = c(200, 470))
hist(nchar(aTNF_ASV), breaks = 20,  xlim = c(200, 470))
library("Biostrings")
pA <- stringDist(c(TRIM_ASV, aTNF_ASV))
saveRDS(pA, "stringDist_ASV.RDS")
