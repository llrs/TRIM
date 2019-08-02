library("WGCNA")
library("RGCCA")
library("integration")
library("org.Hs.eg.db")
folder1 <- "../intestinal_16S_RNAseq_metadb"
RNAseq <- readRDS(file.path(folder1, "expr.RDS"))
otus <- readRDS(file.path(folder1, "otus_table.RDS"))

# model 3
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))


expr_scaled <- scale(RNAseq)
texpr_scaled <- t(expr_scaled)
weight_expr <- t(RNAseq) * drop(model2.2$astar$RNAseq)

otus_f <- otus[rownames(model2.2$a$`16S`), ]
weight_otus <- t((otus_f)) * drop(model2.2$astar$`16S`)

dna_orig <- model2.2$a$`16S`[, 1]
rna_orig <- model2.2$a$RNAseq[, 1]

dna_orig_names <- names(dna_orig)[dna_orig != 0]
rna_orig_names <- names(rna_orig)[rna_orig != 0]


wes <- weight_expr[, rna_orig_names]
wos <- weight_otus[, dna_orig_names]
n_unique <- function(x){length(unique(x))}
sd_rna <- apply(wes, 2, n_unique)
sd_dna <- apply(wos, 2, n_unique)
wes <- wes[, sd_rna != 1]
wos <- wos[, sd_dna != 1]


# Test if there are two groups (one for each sign)
pseudo_pvalue <- function(x) {
  cl <- kmeans(x, centers = 2, iter.max = 100)
  1 - cl$betweenss/cl$totss
}
p_rna <- apply(wes, 2, pseudo_pvalue)
p_dna <- apply(wos, 2, pseudo_pvalue)

# We could use it to remove some correlations as there are two groupsonly it will be a perfect correlation
wes <- wes[, !p_rna == 0]
wos <- wos[, !p_dna == 0]

# Correlation expression (wes) and otus (wos)
cors <- WGCNA::cor(wes, wos, method = "pearson", use = "pairwise.complete.obs")

n_cor <- 158
cor_sign <- cor_sign(n_cor)
df <- as.data.frame(which(abs(cors) > cor_sign, arr.ind = TRUE))
df$cols <- colnames(cors)[df$col]

tax <- readRDS("../intestinal_16S_RNAseq_metadb/otus_tax.RDS")
df$Microorganism <- apply(tax[df$cols, c("Genus", "Species")], 1, paste, collapse = " ")
df$Microorganism <- gsub("NA", "", df$Microorganism)
df$Microorganism <- gsub("^ $", "", df$Microorganism)
df$Genes <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(df)), 
                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
df$cor <- apply(df[, c("row", "col")], 1, function(x, cors){cors[x[1], x[2]]}, cors = cors)
df$pvalue <- sapply(df$cor, pvalue, n = n_cor)
df <- df[abs(order(df$cor)), ]
df$Microorganism <- trimws(df$Microorganism)
# write.csv(df[, c(-1, -2, -3)], file = "correlations_bootstrapping.csv", 
          # row.names = FALSE)

# Whatch out the dispersion of the weights
sds2 <- apply(weight_expr, 2, sd)
means2 <- colMeans(weight_expr)
plot(sds2, means2)
hist(sds2, breaks = 1000)

sds1 <- apply(weight_otus, 2, sd)
means1 <- colMeans(weight_otus)
plot(sds1, means1)
hist(sds1, breaks = 50)

# Probably it needs a filter on there too, 
df_sds2 <- sds2[colnames(weight_expr)[df$row]]
if (sum(!is.na(names(df_sds2))) < nrow(df)) {
  stop("incorrect names")
}
df_sds1 <- sds1[df$cols]
head(df[df_sds2 > 0.005 & df_sds1 > 0.005, ])
plot(weight_expr[, "ENSG00000155265.10"], weight_otus[, "OTU_117"])

