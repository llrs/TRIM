
jaccard <- read.delim("partek_Jaccard_index_dissimilarity_matrix.txt", check.names = FALSE)
cmd <- cmdscale(jaccard[, 2:ncol(jaccard)])

summ <- read.delim("Partek_Microbiome2_Kraken_Classified_species.txt")
unclassified <-  summ[, 3]
pur <- summ[, 3]/(summ[, 2] + summ[, 3])*100

o <- order(pur, decreasing = TRUE)
pur_ord <- pur[o]
cor(pur[o], cmd[, 1][o], method = "spearman")
barplot(pur_ord, main = "Percentage of unclassified reads")
abline(h = 50, col = "red", lwd = 3)

no_rank <- summ[, 12]
o2 <- order(no_rank, decreasing = TRUE)
barplot(no_rank[o2])
plot(no_rank, summ[, 3])
cor.test(no_rank, summ[, 3])
