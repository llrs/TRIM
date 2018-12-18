library("integration")
library("ggplot2")
library("ggraph")
library("tidygraph")
library("ggnetwork")
library("data.table")

files <- list.files(pattern = "toppfun_.*.txt.tsv")
# colC <- rep("character", 11)
# f1 <- read.table(files[[1]], sep = "\t", header = FALSE, quote = "", 
#                  stringsAsFactors = FALSE, dec = ",")
# 

paths <- sapply(files, function(x){
  f1 <- read.table(x, sep = "\t", header = TRUE, quote = "", 
                   stringsAsFactors = FALSE, dec = ",")
  f1[, 3]
})

pvalues <- sapply(files, function(x){
  f1 <- read.table(x, sep = "\t", header = TRUE, quote = "", 
                   stringsAsFactors = FALSE, dec = ",")
  f1[, 5]
})

capital <- function(s){
  paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
        sep = "")
}

pvalues <- as.numeric(unlist(pvalues, use.names = FALSE))
pvalues[pvalues > 1] <- 1
names(paths) <- gsub("toppfun_(.*)\\.txt.tsv", "\\1", names(paths))
names(paths) <- gsub("Methylo", "Methylobacterium", names(paths),
                     ignore.case = TRUE)

df <- cbind.data.frame("from" = rep(capital(names(paths)), lengths(paths)), 
      "to" = unlist(paths, use.names = FALSE),
      "pval" = pvalues)
graph <- as_tbl_graph(df)

tp <- table(df[, 2])
tm <- table(df[, 1])


graph %<>% activate(nodes) %>% 
  mutate(Type = ifelse(name %in% capital(names(paths)), "Microorganism", "Pathway"),
         Connections = ifelse(name %in% names(tp), tp[name], tm[name]))

ggraph(graph, layout = "lgl") +
  # scale_edge_colour_gradient2() +
  geom_edge_link(aes(color = -log10(pval))) +
  geom_node_point(aes(filter = Type != "Microorganism" & Connections >= 2, size = Connections)) +
  # geom_node_point(aes(filter = Type == "Microorganism")) +
  geom_node_label(aes(filter = Type == "Microorganism", label = name),
                  size = 4, repel = TRUE) +
  scale_size(range = c(0.5, 3)) +
  scale_edge_colour_continuous(low = "#56B1F7", high = "#132B43") +
  theme_blank()

a <- split(df$from, df$to)
suba <- a[names(tp)[tp >= 3]]
suba <- suba[order(lengths(suba), decreasing = TRUE)]
seq.max <- seq_len(max(lengths(suba)))
mat <- t(sapply(suba, "[", i = seq.max))

fwrite(as.data.frame(mat), file = "paths2microorganisms.tsv", na = "")
