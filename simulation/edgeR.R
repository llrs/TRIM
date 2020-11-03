library("edgeR")

thinCounts
raw <- readRDS("intestinal_16S_RNAseq_integration/raw_otus.RDS")

m <- matrix(0, nrow = nrow(raw), ncol = ncol(raw))
nrow <- nrow(raw)
for (col in seq_len(ncol(raw))) {
  m[, col] <- rpois(nrow, mean(raw[, col]))
}
boxplot(colMeans(raw), colMeans(m))
boxplot(rowMeans(raw), rowMeans(m)) # Not similar!!


# thinCounts(m) # Freezes the computer
