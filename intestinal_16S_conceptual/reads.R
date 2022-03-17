# To compare the amount of read per site

library("ggplot2")
library("patchwork")

folder <- "../intestinal_16S"
files <- list.files(pattern = "KrakenSummary", path = folder, full.names = TRUE)

report <- lapply(files, read.delim)
o <- lapply(report, function(x){
  columns <- x[, 1:3]
  cbind.data.frame(columns, reads = rowSums(x[, 2:3]))

})
theme_set(theme_bw())
theme_update(axis.text.x = element_blank(),
             panel.grid.major.x = element_blank(),
             axis.ticks.x = element_blank())

results <- sapply(o, function(x){max(x[, 4])})
max_value <- max(results[1:3])
Macrogen <- ggplot(o[[1]]) +
  geom_col(aes(Sample.name, reads)) +
  ylim(0, max_value) +
  labs(title = "Macrogen", x = "Samples", y = "")
Michigan <- ggplot(o[[2]]) +
  geom_col(aes(Sample.name, reads)) +
  ylim(0, max_value) +
  labs(title = "Michigan", x = "", y = "")
Munich <- ggplot(o[[3]]) +
  geom_col(aes(Sample.name, reads)) +
  ylim(0, max_value) +
  labs(title = "Munich", x = "", y = "Reads")
Munich+Macrogen+Michigan
ggsave("Figures/reads.png")


r <- read.delim("reads.txt", header = FALSE, stringsAsFactors = FALSE)
f <- r[c(TRUE, FALSE), ]
rr <- r[rev(c(TRUE, FALSE)), ]
df <- data.frame(files = f, reads = as.numeric(rr))
df <- df[df$reads != 0, ]
pdf("Figures/read_sequences.pdf")
barplot(sort(df$reads), ylim = c(0, 10^5))
barplot(log10(sort(df$reads)), ylim = c(0, 5))
dev.off()
