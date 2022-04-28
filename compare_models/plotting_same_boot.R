library("ggplot2")
library("gplots")
library("dplyr")
library("org.Hs.eg.db")
library("integration")
library("patchwork")

if (basename(getwd()) == "TRIM") {
  folder <- "compare_models/"
  folder0 <- "intestinal_16S_RNAseq_integration"
  folder1 <- "intestinal_16S_RNAseq_metadb"
} else {
  folder <- "."
  folder0 <- "../intestinal_16S_RNAseq_integration"
  folder1 <- "../intestinal_16S_RNAseq_metadb"
}
b0 <- readRDS(file.path(folder, "boot0.RDS"))
b1.2 <- readRDS(file.path(folder, "boot1.2.RDS"))
b2.2 <- readRDS(file.path(folder, "boot2.2.RDS"))

# Models ###
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))

# model 1 without interaction
# model 2 without interaction
model1.2 <- readRDS(file.path(folder1, "model2_best.RDS"))
# model 3
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))

RNAseq <- readRDS(file.path(folder1, "expr.RDS"))


theme_set(theme_bw())
theme_update(strip.background = element_blank(), 
             panel.grid.minor = element_blank())

# AVE ####
AVE0 <- t(sapply(b0, "[[", "AVE"))
AVE1.2 <- t(sapply(b1.2, "[[", "AVE"))
AVE2.2 <- t(sapply(b2.2, "[[", "AVE"))
b <- rbind.data.frame(cbind.data.frame(AVE0, model = "0"), 
                      cbind.data.frame(AVE1.2, model = "1.2"),
                      cbind.data.frame(AVE2.2, model = "2.2"))
b$index <- rep(seq_len(1000), 3)
# This index doesn't converge
b <- b[!(b$inner == 0 & b$outer == 0), ]
ggplot(b) + 
  geom_density(aes(inner, group = model, fill = model), alpha = 0.5)
ggplot(b) + 
  geom_density(aes(outer, group = model, fill = model), alpha = 0.5)

AVE_names <- c("AVE_inner", "AVE_outer")
b$model <- as.character(b$model)
p <- ggplot(b) +
  geom_point(aes(inner, outer, col = model, shape = model), alpha = 0.5) +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model2.2$AVE[AVE_names])[1, , drop = FALSE], 
             fill = "blue", col = "black", shape = 21) +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model1.2$AVE[AVE_names])[1, , drop = FALSE], 
             fill = "green", col = "black", shape = 21) +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model0$AVE[AVE_names])[1, , drop = FALSE], 
             fill = "red", col = "black", shape = 21) +
  stat_ellipse(aes(inner, outer, col = model)) +
  labs(title = "AVE in bootstraps in HSCT dataset", x = "Inner AVE", y = "Outer AVE", col = "Model", shape = "Model") +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.01))

ggsave("Figures/hsct_bootstrap_models.png", plot = p, dpi = 300, units = "mm")
ggsave("Figures/Figure_6.png", plot = p, width = 170, dpi = 300, units = "mm")
ggsave("Figures/Figure_6_half.png", plot = p, width = 85, dpi = 300, units = "mm")

b %>% 
  dplyr::select(-index) %>% 
  tidyr::gather(inner, outer, -model) %>% 
  dplyr::rename(AVE = inner, value = outer) %>% 
  group_by(model, AVE) %>% 
  summarise(mean(value), sd(value)) %>% 
  arrange(AVE) %>% 
  write.csv(file = "dispersion_models.csv", row.names = FALSE)

index <- readRDS("compare_models/index_locale.RDS")
meta_r <- readRDS("intestinal_16S_RNAseq_metadb/meta.RDS")
meta <- meta_r

i <- sapply(index, function(x)x)
t_i <- sort(table(i)) - 1000
barplot(t_i, main = "Deviation from 1000", 
        col = as.factor(meta$Patient_ID[as.numeric(names(t_i))]))
barplot(t_i, main = "Deviation from 1000", 
        col = as.factor(meta$IBD[as.numeric(names(t_i))]))

df <- as.data.frame(table(i) - 1000)
meta2 <- cbind(df, meta)

ggplot(meta2) + 
  geom_col(aes(forcats::fct_reorder(i, Freq), Freq, col = AGE_SAMPLE, fill = AGE_SAMPLE)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.border = element_blank()) +
  labs(x = "Sample", y = "Deviation from 1000", 
       title = "Samples repeated in bootstrapping",
       fill = "Age", col = "Age") +
  scale_color_viridis_c(aesthetics = c("fill", "colour"), direction = -1)
ggplot(meta2) + 
  geom_col(aes(forcats::fct_reorder(i, Freq), Freq, col = IBD, fill = IBD)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.border = element_blank()) +
  labs(x = "Sample", y = "Deviation from uniform distribution", 
       title = "Samples repeated in bootstrapping",
       col = "Disease", fill = "Disease")
ggplot(meta2) + 
  geom_col(aes(forcats::fct_reorder(i, Freq), Freq, col = SEX, fill = SEX)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.border = element_blank()) +
  labs(x = "Sample", y = "Deviation from uniform distribution", 
       title = "Samples repeated in bootstrapping",
       col = "Sex", fill = "Sex")


df <- sapply(index, function(i, x) {
  
  c(Controls = sum(x$IBD[i] != "CD")/158,
    Ileum = sum(x$Exact_location[i] == "ILEUM", na.rm = TRUE)/158,
    Smoker = sum(x$Tobacco[i] == "No", na.rm = TRUE)/158,
    Age = mean(x$AGE_SAMPLE[i]),
    Female = sum(x$SEX == "FEMALE")/158
  )
}, x = meta_r)

df <- as.data.frame(t(df))
df2 <- data.frame(Controls = sum(meta_r$IBD != "CD")/158,
                  Ileum = sum(meta_r$Exact_location == "ILEUM", na.rm = TRUE)/158,
                  Smoker = sum(meta_r$Tobacco == "No", na.rm = TRUE)/158,
                  Age = mean(meta_r$AGE_SAMPLE),
                  Female = sum(meta_r$SEX == "FEMALE")/158)

CC <- ggplot(df) +
  geom_count(aes(Controls, Ileum), col = "grey") +
  geom_point(aes(Controls, Ileum), col = "black", data = df2) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  # scale_x_continuous(labels = scales::percent) +
  labs(
    # x = "Controls (%)",  
    y = "Ileum (%)",
    n = "Bootstraps",
    title = "Bootstrap samples",
    subtitle = paste0(length(index), " resamples of all the samples"))

CA <- df %>% 
  ggplot(aes(Controls, Age)) +
  geom_point(col = "grey") +
  geom_smooth( col = "darkgrey") +
  geom_point(col = "black", data = df2) +
  theme_minimal() +
  # scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Controls (%)",  
       y = "Age (Mean)")
plot_boot <- CC/CA
plot_boot
ggsave("Figures/Figure9.png", dpi = 300, width = 170, units = "mm")

# microbiome ####
micro_fun <- function(x){x$STAB$`16S`}
micro2.2 <- lapply(b2.2, micro_fun)
micro2.2 <- micro2.2

dna <- matrix(0, nrow = nrow(model0$a$`16S`), ncol = 1000)
rownames(dna) <- rownames(model0$a$`16S`)
i <- 1
for (mod in micro2.2) {
  dna[, i][names(mod)] <- mod
  i <- i + 1
}

# transcriptome ####
rna_fun <- function(x){x$STAB$RNAseq}
rna2.2 <- lapply(b2.2, rna_fun)
rna2.2 <- rna2.2

rna <- matrix(0, nrow = nrow(RNAseq), ncol = 1000)
rownames(rna) <- rownames(RNAseq)
i <- 1
for (mod in rna2.2) {
  rna[, i][names(mod)] <- mod
  i <- i + 1
}

# Correlations ####
cors <- WGCNA::cor(t(rna), t(dna), 
                   method = "spearman", 
                   use = "pairwise.complete.obs")
saveRDS(cors, "cors_bootstrap.RDS")
dna_orig <- model2.2$a$`16S`[, 1]
rna_orig <- model2.2$a$RNAseq[, 1]

dna_orig_names <- names(dna_orig)[dna_orig != 0]
rna_orig_names <- names(rna_orig)[rna_orig != 0]

sub_cors2 <- WGCNA::cor(t(rna[rna_orig_names, ]), 
                        t(dna[dna_orig_names, ]), 
                        method = "spearman", 
                        use = "pairwise.complete.obs")
sub_cors <- subcors2
sub_cors <- cors[rownames(cors) %in% rna_orig_names, 
                 colnames(cors) %in% dna_orig_names]

heatmap.2(sub_cors, scale = "none", labCol = FALSE,
          labRow = FALSE, xlab = "Microorganisms (OTUs)", ylab = "Transcripts",
          main = "Correlations of bootstrapping model 2.2", margins = c(2, 2), 
          col = scico::scico(50, palette = 'roma'), trace = "none")

# Test if there are two groups (one for each sign)
pseudo_pvalue <- function(x) {
  cl <- kmeans(x, centers = 2, iter.max = 100)
  1 - cl$betweenss/cl$totss
  }
p_rna <- apply(rna, 1, pseudo_pvalue)
p_dna <- apply(dna, 1, pseudo_pvalue)
hist(p_rna, breaks = 1000)
hist(p_rna[rna_orig_names], breaks = 1000)
non_grouped <- names(p_rna)[p_rna < 0.4 & p_rna > 0.3]

n_cor <- 1000
cor_sign <- cor_sign(n_cor)
df <- as.data.frame(which(sub_cors > cor_sign, arr.ind = TRUE))
df$cols <- colnames(sub_cors)[df$col]

tax <- readRDS("../intestinal_16S_RNAseq_metadb/otus_tax.RDS")
df$Microorganism <- apply(tax[df$cols, c("Genus", "Species")], 1, paste, collapse = " ")
df$Microorganism <- gsub("NA", "", df$Microorganism)
df$Microorganism <- gsub("^ $", "", df$Microorganism)
df$Genes <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(df)), 
                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
df$cor <- apply(df[, c("row", "col")], 1, function(x, cors){cors[x[1], x[2]]}, cors = sub_cors)
df$pvalue <- sapply(df$cor, pvalue, n = n_cor)
df2 <- df[rownames(df) %in% non_grouped, ]
write.csv(df[, c(-1, -2, -3)], file = "correlations_bootstrapping.csv", 
          row.names = FALSE)
