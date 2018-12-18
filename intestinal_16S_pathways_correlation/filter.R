# Load ####
cd <- setwd("..")

intestinal <- "intestinal_16S"

file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
rna <- "intestinal_RNAseq"
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)


setwd(cd)

# Libraries
library("integration")
library("org.Hs.eg.db")
library("reactome.db")

today <- format(Sys.time(), "%Y%m%d")

gsv <- readRDS("gsv.RDS")
genus_i <- readRDS("genus.RDS")

meta_r <- meta_r_norm(meta_r)

all_s <- readRDS("correlations_all.RDS")
pall_s <- readRDS("padj_all.RDS")

filter_sexual <- function(expr) {
  source("../genes_XY.R")
  sexual_related <- gsub("(.+)\\..*", "\\1", rownames(expr)) %in% 
  c(#bmX$ensembl_gene_id, 
    bmY$ensembl_gene_id)
  expr[!sexual_related, ]
}

# Output ####
#  The numbers come from the number of samples in each correlation
threshold <- 0.05
all_samples_ensembl <- relevant("../intestinal_16S_pathways_metadb/sgcca.RDS", all_s, pall_s, threshold)
all_samples_symbol <- all_samples_ensembl
all_samples_symbol$Gene <- mapIds(reactome.db, keys = all_samples_symbol$Gene, 
                                  keytype = "PATHID", column = "PATHNAME")
all_samples_symbol$Gene <- gsub("Homo sapiens: ", "", all_samples_symbol$Gene)
all_samples_symbol <- all_samples_symbol[!is.na(all_samples_symbol$Gene), ]
all_samples_symbol <- all_samples_symbol[!duplicated(all_samples_symbol), ]
write.csv(all_samples_symbol, file = "correlation_all.csv", row.names = FALSE, na = "")

# 
# ibd_ensembl <- relevant("IBD.RData", disease, pdisease, threshold)
# ibd_symbol <- ibd_ensembl
# ibd_symbol$Gene <- gsub("(.*)\\..*", "\\1", ibd_symbol$Gene)
# ibd_symbol$Gene <- mapIds(org.Hs.eg.db, keys = ibd_symbol$Gene, 
#                                   keytype = "ENSEMBL", column = "SYMBOL")
# ibd_symbol <- ibd_symbol[!is.na(ibd_symbol$Gene), ]
# write.csv(ibd_symbol, file = "correlation_CD.csv", row.names = FALSE, na = "")
# 
# contr <- relevant("Controls.RData", controls, pcontrols, threshold)
# contr$Gene <- gsub("(.*)\\..*", "\\1", contr$Gene)
# contr$Gene <- mapIds(reactome.db, keys = contr$Gene, keytype = "PATHID", column = "PATHNAME")
# write.csv(contr, file = "correlation_Controls.csv", row.names = FALSE, na = "")

plot_single_cor <- function(x, gene, y, microorganism, colr, case, cor_val) {
  # genes <- gsub("(.*)\\..*", "\\1", gene)
  symbol <- tryCatch({mapIds(
    reactome.db, key = gene, keytype = "PATHID",
    column = "PATHNAME"
  )},  error = function(e){NA}) 

    if (is.na(symbol)) {
    return(NA)
  }
  
  x_s <- x[gene, ]
  x_s[x_s == 0] <- NA
  
  y_s <- y[microorganism, ]
  y_s[y_s == 0] <- NA
  pch_start <- 15
  
  # main <- cor(x_s, y_s, method = "spearman", use = "pairwise.complete.obs")
  plot(x_s, y_s, xlab = paste(symbol, collapse = " "), ylab = microorganism, main = cor_val, col = colr, 
       pch = pch_start + as.numeric(case))
  legend("bottomleft", fill = as.factor(levels(colr)), legend = levels(colr))
  legend("topright", pch = pch_start + as.numeric(as.factor(levels(case))), 
         legend = levels(case))
  # legend("topleft", lty = "solid", col = as.factor(levels(bg)), 
         # legend = levels(bg), lwd = 2)
}

pdf("correlations_plot_all.pdf")
o <- apply(all_samples_ensembl, 1, function(x) {
  micro <- x[[1]]
  gene <- x[[2]]
  cor_val <- x[[3]]
  plot_single_cor(gsv, gene, genus_i, micro, as.factor(meta_r$Time), 
                  as.factor(meta_r$Involved_Healthy), 
                  cor_val = cor_val)
})
dev.off()

pdf("correlations_plot_ibd.pdf")
o <- apply(ibd_ensembl, 1, function(x) {
  micro <- x[[1]]
  gene <- x[[2]]
  cor_val <- x[[3]]
  keep <- meta_r$IBD != "CONTROL"
  plot_single_cor(gsv[, keep], gene, genus_i[, keep], micro, as.factor(meta_r$Exact_location[keep]), 
                  as.factor(meta_r$SEX[keep]), cor_val = cor_val)
})
dev.off()

## Network representation
# FIXME: Loaded but now it is not recognized as a network
# https://www.biostars.org/p/80498/
# bact <- all_samples_ensembl$Gene[all_samples_ensembl$Microorganism == "Bacteroides"]
# M <- disease[, bact]
# edges <- NULL
# for (i in 1:nrow(M)) {
#   for (j in 1:ncol(M)) {
#     edges <- rbind(edges, c(rownames(M)[i], colnames(M)[j], M[i,j]))
#   }
# }
# 
# colnames(edges) <- c('node1', 'node2', 'value')
# write.table(edges, 'edges.txt', row.names=FALSE, quote=FALSE, sep='\t')
# 
# library("ggplot2")
# library("ggraph")
# library("tidygraph")
# library("ggnetwork")
# tab_names <- table(all_samples_symbol$Microorganism) > 60
# org <- names(tab_names[tab_names])
# org <- org[!org %in% "Peptostreptococcus"]
# 
# graph_data <- all_samples_symbol[all_samples_symbol$Microorganism %in% org, ]
# 
# tg <- table(graph_data$Gene[graph_data$Microorganism %in% org])
# tm <- table(graph_data$Microorganism[graph_data$Microorganism %in% org])
# t_number <- c(tg, tm)
# 
# colnames(graph_data) <- c("from", "to", "Correlation", "pvalue")
# tidyg <- as_tbl_graph(graph_data, directed = FALSE)
# 
# # Add info about the nodes
# tidyg2 <- tidyg %>% 
#   activate(nodes) %>% 
#   mutate(Connections = t_number[name], 
#          Type = ifelse(name %in% org, "Microorganism", "Gene"),
#          Sth = c(""))
# 
# p <- ggraph(tidyg2, layout = "lgl") +
#   scale_edge_colour_gradient2() +
#   geom_edge_link(aes(colour = Correlation)) +
#   geom_node_point(aes(filter = Type != "Microorganism" & Connections >= 2, size = Connections)) +
#   # geom_node_point(aes(filter = Type == "Microorganism", size = Connections)) +
#   geom_node_label(aes(filter = Type == "Microorganism", label = name), 
#                   size = 4, repel = TRUE) +
#   scale_size(range = c(0.5, 2)) +
#   theme_blank()
# ggsave("Figures/connections.png")
  