library("ggplot2")
library("patchwork")

folder <- "../intestinal_16S"
files <- list.files(pattern="KrakenSummary", path = folder, full.names = TRUE)

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
  labs(title = "Macrogen", x = "", y = "Reads")
Michigan <- ggplot(o[[2]]) +
  geom_col(aes(Sample.name, reads)) +
  ylim(0, max_value) +
  labs(title = "Michigan", x = "Samples", y = "")
Munich <- ggplot(o[[3]]) +
  geom_col(aes(Sample.name, reads)) +
  ylim(0, max_value) +
  labs(title = "Munich", x = "", y = "")

Macrogen + Michigan + Munich
