library("ggplot2")

folder <- "." # intestinal_16S_RNAseq_integration
b0 <- readRDS(file.path(folder, "boot0.RDS"))
b1.2 <- readRDS(file.path(folder, "boot1.2.RDS"))
b2.2 <- readRDS(file.path(folder, "boot2.2.RDS"))

# Models ###
folder0 <- "../intestinal_16S_RNAseq_integration"
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))

# model 1 without interaction
folder1 <- "../intestinal_16S_RNAseq_metadb"
# model 2 without interaction
model1.2 <- readRDS(file.path(folder1, "model2_best.RDS"))
# model 3
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))

RNAseq <- readRDS(file.path(folder1, "expr.RDS"))

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
ggplot(b) +
  geom_point(aes(inner, outer, col = model), alpha = 0.5) +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model2.2$AVE[AVE_names]), 
             col = "blue") +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model1.2$AVE[AVE_names])[1, , drop = FALSE], 
             col = "green") +
  geom_point(aes(AVE_inner, AVE_outer), 
             data = as.data.frame(model0$AVE[AVE_names])[1, , drop = FALSE], 
             col = "red") +
  labs(title = "AVE in bootstraps")

i <- sapply(index, function(x)x)
t_i <- sort(table(i)) - 1000
barplot(t_i, main = "Deviation from 1000", 
        col = as.factor(meta$Patient_ID[as.numeric(names(t_i))]))
barplot(t_i, main = "Deviation from 1000", 
        col = as.factor(meta$IBD[as.numeric(names(t_i))]))

df <- as.data.frame(table(i)-1000)
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

meta_r <- readRDS("../intestinal_16S_RNAseq_metadb//meta.RDS")
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
        axis.title.x = element_blank()) +
  # scale_x_continuous(labels = scales::percent) +
  labs(
    # x = "Controls (%)",  
       y = "Ileum (%)",
       title = "Distribution of the bootstrapping samples",
       subtitle = paste0(length(index), " resamples of 158 samples"))

CA <- ggplot(df) +
  geom_point(aes(Controls, Age), col = "grey") +
  geom_smooth(aes(Controls, Age), col = "darkgrey") +
  geom_point(aes(Controls, Age), col = "black", data = df2) +
  theme_minimal() +
  # scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Controls (%)",  
       y = "Age (Mean)",
       caption = "HSCT cohort")
CC/CA

# microbiome ####
micro_fun <- function(x){x$STAB$`16S`}
micro2.2 <- lapply(b2.2, micro_fun)
micro2.2 <- micro2.2[-230]

dna <- matrix(0, nrow = nrow(model0$a$`16S`), ncol = 999)
rownames(dna) <- rownames(model0$a$`16S`)
i <- 1
for (mod in micro2.2) {
  dna[, i][names(mod)] <- mod
  i <- i + 1
}

# transcriptome ####
rna_fun <- function(x){x$STAB$RNAseq}
rna2.2 <- lapply(b2.2, rna_fun)
rna2.2 <- rna2.2[-230]

rna <- matrix(0, nrow = nrow(RNAseq), ncol = 999)
rownames(rna) <- rownames(RNAseq)
i <- 1
for (mod in rna2.2) {
  rna[, i][names(mod)] <- mod
  i <- i + 1
}

cors <- WGCNA::cor(t(rna), t(dna), 
                   method = "spearman", 
                   use = "pairwise.complete.obs")

cors <- cor(t(rna), t(dna), use = "pairwise.complete.obs")
