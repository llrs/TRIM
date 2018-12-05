library("gghighlight")
library("org.Hs.eg.db")
library("ggplot2")
library("ggforce")
library("RGCCA")
library("Matrix")
library("tidyr")
library("dplyr")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")

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

# Some groups are like clusters
plot(Ile_CD_vs_Co_Ctr, Co_CD_vs_Ile_Ctr)
abline(a = 0, b = 1, col = "red")

hist(Ile_CD_vs_Co_Ctr-Co_CD_vs_Ile_Ctr)
# Normal distribution

d <- cbind.data.frame(ICD_CCt = Ile_CD_vs_Co_Ctr,
                      CCD_ICt = Co_CD_vs_Ile_Ctr,
                      CCD_CCt = Co_CD_vs_Co_Ctr,
                      ICD_ICt = Ile_CD_vs_Ile_Ctr,
                      ICD_CCD = Co_CD_vs_Ile_CD,
                      CCt_ICt = Co_Ctr_vs_Ile_Ctr)
saveRDS(d, "distance_centroids.RDS")

ggplot(d) +
  geom_density2d(aes(ICD_CCt, CCD_ICt)) +
  geom_abline(slope = 1, intercept = 0)
ggplot(d) +
  geom_point(aes(ICD_CCt, CCD_ICt), alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, col = "red")
d2 <- d[, 1] - d[, 2]
d4 <- d[, 3] - d[, 4]

# Ok we believe that those around distance d2 ~= 0 are the nice ones
# Also the ones with more distance between themselfs could be interesting
# They follow like a pareto distribution:
hist(rowSums(d))

# Now retrieve the names!!
picky <- pairs[abs(d2) < 0.002] # To be picky and "only" ~800 pairs
pic <- simplify2array(picky)

taxa <- otus_tax_i[pic[1, ], 6:7]
taxa[is.na(taxa)] <- ""
taxes <- apply(taxa, 1, paste0, collapse = " ")

g <- trimVer(pic[2, ])
g2 <- mapIds(org.Hs.eg.db, keys = g, keytype = "ENSEMBL", column = "SYMBOL")

paired <- cbind.data.frame("Micro" = taxes, "Gene" = g2)
saveRDS(paired, "gene_microorganism.RDS")

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
ggplot(denrich) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = Micro), alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  # gghighlight(ICD_CCt > 10) +
  guides(color = FALSE)

droplets <- denrich %>%
  filter(ICD_CCt > 7.5 & CCD_ICt >7.5) %>% 
  group_by(Gene) %>% 
  summarise(ICD_CCt_max = max(ICD_CCt),
            ICD_CCt_min = min(ICD_CCt),
            CCD_ICt_max = max(CCD_ICt),
            CCD_ICt_min = min(CCD_ICt)
            )
droplets2 <- denrich %>% 
  filter(ICD_CCt > 7.5 | CCD_ICt > 7.5)
write.csv(droplets2, file = "distancies_meteorits.csv")

a <- denrich %>% 
  group_by(Gene) %>% 
  count(n_distinct(Ensembl)) %>% 
  ungroup()

keep <- denrich$ICD_CCt > 7.5 | denrich$CCD_ICt > 7.5
keep_f <- function(x){
  denrich$ICD_CCt > x | denrich$CCD_ICt > x
}
sapply(1:13, function(x){length(unique(denrich$Gene[keep_f(x)]))})
ggplot(droplevels(denrich)) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = "grey"), alpha = 0.2) +
  geom_point(data = denrich[keep, ], aes(ICD_CCt, CCD_ICt, color = Gene)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  guides(color = FALSE)
keep2 <- denrich$ICD_CCt > 10.5 & denrich$CCD_ICt < 5
ggplot(droplevels(denrich)) +
  geom_point(aes(ICD_CCt, CCD_ICt, color = "grey"), alpha = 0.2) +
  geom_point(data = denrich[keep2, ], aes(ICD_CCt, CCD_ICt, color = Gene)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  guides(color = FALSE)

rank <- rowSums(denrich[, 1:6])
library("fgsea")
names(rank) <- as.character(seq_along(rank))
nam <- split(seq_along(rank), denrich$Micro)
nam <- nam[-which(names(nam) == " ")] # Remove the empty one
microEnrich <- fgsea(nam, stats = rank, nperm = 10000)
microEnrich %>% 
  filter(padj < 0.05) %>% 
  arrange(padj, desc(abs(NES))) %>% 
  select(-leadingEdge, -size) %>% 
  as_tibble()

nam2 <- split(seq_along(rank), denrich$Gene)
geneEnrich <- fgsea(nam2, stats = rank, nperm = 10000)
geneEnrich %>% 
  filter(padj < 0.05) %>% 
  arrange(padj, desc(abs(NES))) %>% 
  select(-leadingEdge, -size) %>% 
  as_tibble()
si <- geneEnrich %>% 
  filter(padj < 0.05)
hist(si$NES)
