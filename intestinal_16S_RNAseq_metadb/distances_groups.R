library("patchwork")
library("gghighlight")
library("org.Hs.eg.db")
library("ggplot2")
library("ggforce")
library("RGCCA")
library("Matrix")
library("tidyr")
library("dplyr")
library("integration")
library("fgsea")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("otus_table.RDS")
otus_tax_i <- readRDS("otus_tax.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS("meta.RDS")


meta_r$Location <- ifelse(meta_r$Exact_location != "ILEUM", "COLON", "ILEUM")
levels(as.factor(meta_r$Location))
levels(as.factor(meta_r$IBD))
keep <- allComb(meta_r, c("Location", "IBD"))

model <- readRDS("model3_best.RDS")

otus3 <- rownames(model3$a[[2]])[model3$a[[2]][, 1] != 0 ]
otus2 <- rownames(model2$a[[2]])[model2$a[[2]][, 1] != 0 ]
o <- intersect(otus3, otus2)

relevantOTUs <- rownames(model$a[[2]])[model$a[[2]][, 1] != 0]
relevantGenes <- rownames(model$a[[1]])[model$a[[1]][, 1] != 0]
otus <- otus_table_i[relevantOTUs, ]
expr_f <- expr[relevantGenes, ]

# Calculate the center of each group
out <- lapply(as.data.frame(keep), function(x) {
  Y <- rowMeans(otus[, x], na.rm = TRUE)
  X <- rowMeans(expr_f[, x], na.rm = TRUE)
  list("Y" = Y, "X" = X)
})

m <- length(relevantOTUs)*length(relevantGenes)
apoint <- function(x, m) {
  out_points <- vector("list", length = m)
  n <- 1
  for (o in x$Y) {
    for (i in x$X) {
      out_points[[n]] <- c(o, i)
      n <- n+1
    }
  }
  out_points
}
# Join the names of the pairs genes otus
aname <- function(otus, expr) {
  m <- nrow(otus) * nrow(expr)
  out_pairs <- vector("list", length = m)
  n <- 1
  for (o in rownames(otus)) {
    for (i in rownames(expr)) {
      out_pairs[[n]] <- c(o, i)
      n <- n+1
    }
  }
  out_pairs
}
pairs <- aname(otus, expr_f)

points <- lapply(out, apoint, m = m)
p <- simplify2array(points)
redist <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  z <- x-y
  sqrt(sum(z^2))
}

Co_CD_vs_Ile_Ctr <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 1], p[x, 4])
})
Ile_CD_vs_Co_Ctr <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 2], p[x, 3])
})
Co_CD_vs_Ile_CD <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 1], p[x, 2])
})
Co_Ctr_vs_Ile_Ctr <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 3], p[x, 4])
})
Co_CD_vs_Co_Ctr <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 1], p[x, 3])
})
Ile_CD_vs_Ile_Ctr <- sapply(seq_len(nrow(p)), function(x) {
  redist(p[x, 2], p[x, 4])
})

d <- cbind.data.frame(ICD_CCt = Ile_CD_vs_Co_Ctr,
                      CCD_ICt = Co_CD_vs_Ile_Ctr,
                      CCD_CCt = Co_CD_vs_Co_Ctr,
                      ICD_ICt = Ile_CD_vs_Ile_Ctr,
                      ICD_CCD = Co_CD_vs_Ile_CD,
                      CCt_ICt = Co_Ctr_vs_Ile_Ctr)

genes <- vapply(pairs, function(x){
  x[2]
}, character(1L))
tax <- vapply(pairs, function(x){
  x[1]
}, character(1L))

ga <- mapIds(org.Hs.eg.db, keys = trimVer(genes), keytype = "ENSEMBL", column = "SYMBOL")
genus <- otus_tax_i[tax, 6:7]
genus[is.na(genus)] <- ""
genus <- apply(genus, 1, paste0, collapse = " ")
denrich <- cbind(d, Micro = genus, Gene = ga, Ensembl = trimVer(genes))
saveRDS(denrich, "distance_centroids_model3.RDS")

en <- intersect(denrich_model2$Ensembl, denrich_model3$Ensembl)
su <- droplevels(denrich_model3[denrich_model3$Ensembl %in% en, ])
m <- ceiling(max(denrich_model3[, 1:6]))
# Differences between Colon and Ileum
crossed <- ggplot(su) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = Micro), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  # gghighlight(ICD_CCt > 10) +
  guides(color = FALSE) +
  scale_color_viridis_d() +
  xlim(0, m) +
  ylim(0, m) +
  coord_fixed()
# Differences between Controls and Patients
loc <- ggplot(su) +
  geom_point(aes(ICD_ICt, CCD_CCt, color = Micro), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  # gghighlight(ICD_CCt > 10) +
  guides(color = FALSE) +
  xlim(0, m) +
  ylim(0, m) +
  scale_color_viridis_d() +
  coord_fixed()
# Other comparisons
stage <- ggplot(su) +
  geom_point(aes(CCt_ICt, ICD_CCD, color = Micro), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  xlim(0, m) +
  ylim(0, m) +
  # gghighlight(ICD_CCt > 10) +
  guides(color = FALSE) +
  scale_color_viridis_d() +
  coord_fixed()

loc+crossed+stage

o <- sapply(d[, 1:6], ecdf)
plot(o$ICD_CCt)
plot(o$CCD_ICt, add = TRUE, col = "red")
plot(o$CCD_CCt, add = TRUE, col = "blue")
plot(o$ICD_ICt, add = TRUE, col = "green")
plot(o$ICD_CCD, add = TRUE, col = "pink")
plot(o$CCt_ICt, add = TRUE, col = "brown")
legend("bottomright", legend = names(o), 
       fill = c("black", "red", "blue", "green", "pink", "brown"))

isoforms <- su %>% 
  filter(!is.na(Gene)) %>% 
  add_count(Gene) %>% 
  group_by(Gene) %>% 
  summarise(Isoforms = n_distinct(Ensembl)) %>% 
  ungroup()

keep <- denrich$ICD_CCt > 7.5 | denrich$CCD_ICt > 7.5
keep_f <- function(x) {
  denrich$ICD_CCt > x | denrich$CCD_ICt > x
}
sapply(1:13, function(x){length(unique(denrich$Gene[keep_f(x)]))})
ggplot(droplevels(denrich)) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = "grey"), alpha = 0.2) +
  geom_point(data = denrich[keep, ], aes(ICD_CCt, CCD_ICt, color = Ensembl)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  guides(color = FALSE)
keep2 <- denrich$Gene %in% "REG1B"
ggplot(droplevels(denrich)) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = "grey"), alpha = 0.2) +
  geom_point(data = denrich[keep2, ], aes(ICD_CCt, CCD_ICt, color = Gene)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  guides(color = FALSE)
keep3 <- denrich$Gene %in% isoforms$Gene[isoforms$Isoforms != 1 & 
                                           !is.na(isoforms$Gene)]
ggplot(droplevels(denrich[!is.na(denrich$Gene), ])) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = "grey"), alpha = 0.2) +
  geom_point(data = droplevels(denrich[keep3, ]), 
             aes(ICD_CCt, CCD_ICt, color = Ensembl, shape = Gene)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  guides(color = FALSE, shape = FALSE)

rank <- rowSums(denrich_model3[, 1:6])
names(rank) <- as.character(seq_along(rank))
nam <- split(seq_along(rank), denrich$Micro)
nam <- nam[-which(names(nam) == " ")] # Remove the empty one
microEnrich <- fgsea(nam, stats = rank, nperm = 10000)
microEnrich %>% 
  filter(padj < 0.05) %>% 
  arrange(padj, desc(abs(NES))) %>% 
  select(-leadingEdge, -size) %>% 
  View(title = "MicroErnich")
data.table::fwrite(microEnrich, file = "microEnrich_model3_distances_rank.csv")
nam2 <- split(seq_along(rank), denrich$Gene) # The genes without name are discarted
nam3 <- append(nam2, list(intersect = names(rank)[denrich_model3$Gene %in% su$Gene]))
geneEnrich <- fgsea(nam3, stats = rank, nperm = 10000)
geneEnrich %>% 
  filter(padj < 0.05) %>% 
  arrange(padj, desc(abs(NES))) %>% 
  select(-leadingEdge, -size) %>% 
  View(title = "GeneErnich")
si <- geneEnrich %>% 
  filter(padj < 0.05)
hist(si$NES)
data.table::fwrite(geneEnrich, file = "geneEnrich_model3_distances_rank.csv")
