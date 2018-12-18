# Load ####
library("integration")
library("org.Hs.eg.db")
library("reactome.db")
library("ggplot2")
library("ggraph")
library("tidygraph")
library("ggnetwork")
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



today <- format(Sys.time(), "%Y%m%d")

expr <- readRDS("expr.RDS")
meta_r <- droplevels(meta_r[meta_r$`Sample Name_RNA` %in% colnames(expr), ])

# Change each time: species or genus
type <- "genus"
genus_i <- readRDS(paste0(type, ".RDS"))

all_s <- readRDS(paste0("correlations_all_", type,".RDS"))
disease <- readRDS(paste0("correlations_CD_", type,".RDS"))
controls <- readRDS(paste0("correlations_Controls_", type,".RDS"))
colon <- readRDS(paste0("correlations_Colon_", type,".RDS"))
ileum <- readRDS(paste0("correlations_Ileum_", type,".RDS"))

pall_s <- readRDS(paste0("padj_all_", type, ".RDS"))
pdisease <- readRDS(paste0("padj_CD_", type, ".RDS"))
pcontrols <- readRDS(paste0("padj_Controls_", type, ".RDS"))
pcolon <- readRDS(paste0("padj_Colon_", type, ".RDS"))
pileum <- readRDS(paste0("padj_Ileum_", type, ".RDS"))

pre <- "../intestinal_16S_RNAseq_metadb"

filter_sexual <- function(expr) {
  source("../genes_XY.R")
  sexual_related <- gsub("(.+)\\..*", "\\1", rownames(expr)) %in% 
  c(#bmX$ensembl_gene_id, 
    bmY$ensembl_gene_id)
  expr[!sexual_related, ]
}


# Output ####
# #  The numbers come from the number of samples in each correlation
threshold <- 0.05
all_comp <- readSGCCA("../intestinal_16S_RNAseq_integration/sgcca.RDS")
all_samples_ensembl <- relevant(all_comp, all_s, pall_s, threshold)
o2 <- pathsPerMicro(all_samples_ensembl, all_comp)
all_samples_ensembl_cor <- sign_cor(all_s, pall_s, threshold)
write_cor(all_samples_ensembl, file = paste0(today, "_", type, "_correlation_all_model0.csv"))
# write_cor(all_samples_ensembl_cor, file = paste0(today, "_sign_", type, "_correlation_all.csv"))
lapply(names(o2), store_micro, o2 = o2, label = "_model0_all")

colon_samples_ensembl <- relevant(all_comp, colon, pcolon, threshold)
o2 <- pathsPerMicro(colon_samples_ensembl, all_comp)
colon_samples_ensembl_cor <- sign_cor(colon, pcolon, threshold)
write_cor(colon_samples_ensembl, file = paste0(today, "_", type, "_correlation_colon_model0.csv"))

ileum_samples_ensembl <- relevant(all_comp, ileum, pileum, threshold)
o2 <- pathsPerMicro(ileum_samples_ensembl, all_comp)
lapply(names(o2), store_micro, o2 = o2, label = "_model0_ileum")
ileum_samples_ensembl_cor <- sign_cor(ileum, pileum, threshold)
write_cor(ileum_samples_ensembl, file = paste0(today, "_", type, "_correlation_ileum_model0.csv"))


disease_comp <- readSGCCA("../intestinal_16S_RNAseq_integration/IBD.RDS")
ibd_ensembl <- relevant(disease_comp, disease, pdisease, threshold)
ibd_ensembl_cor <- sign_cor(disease, pdisease, threshold)
write_cor(ibd_ensembl, file = paste0(today, "_", type, "_correlation_CD.csv"))
# write_cor(ibd_ensembl_cor, file = paste0(today, "_sign_", type, "_correlation_CD.csv"))


contr_comp <- readSGCCA("../intestinal_16S_RNAseq_integration/Controls.RDS")
contr <- relevant(contr_comp, controls, pcontrols, threshold)
contr_cor <- sign_cor(controls, pcontrols, threshold)
write_cor(contr, file = paste0(today, "_", type, "_correlation_Controls.csv"))
# write_cor(contr_cor, file = paste0(today, "_sign_", type, "_correlation_Controls.csv"))



c_disease <- boots_corr(pdisease, disease_comp)
c_all <- boots_corr(pall_s, all_comp)
c_controls <- boots_corr(pcontrols, contr_comp)


pdf(paste0(today, "_correlations_", type, "_plot_colon.pdf"))
o <- apply(colon_samples_ensembl, 1, function(x) {
  micro <- x[[1]]
  gene <- x[[2]]
  cor_val <- x[[3]]
  plot_single_cor(expr, gene, genus_i, micro, as.factor(meta_r$Time), 
                  as.factor(meta_r$Involved_Healthy), 
                  cor_val = cor_val)
})
dev.off()

pdf(paste0(today, "_correlations_", type, "_plot_ibd.pdf"))
o <- apply(ibd_ensembl, 1, function(x) {
  micro <- x[[1]]
  gene <- x[[2]]
  cor_val <- x[[3]]
  keep <- meta_r$IBD != "CONTROL"
  plot_single_cor(expr[, keep], gene, genus_i[, keep], micro, as.factor(meta_r$Exact_location[keep]), 
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


all_samples_symbol <- ensembl2symbol(all_samples_ensembl)
tab_names <- table(all_samples_symbol$Microorganism) > 10
org <- names(tab_names[tab_names])
org <- org[!org %in% c("Peptostreptococcus", "Dialister")]

graph_data <- all_samples_symbol[all_samples_symbol$Microorganism %in% org, ]

tg <- table(graph_data$Gene[graph_data$Microorganism %in% org])
tm <- table(graph_data$Microorganism[graph_data$Microorganism %in% org])
t_number <- c(tg, tm)

colnames(graph_data) <- c("from", "to", "Correlation", "pvalue")
tidyg <- as_tbl_graph(graph_data, directed = FALSE)

# Add info about the nodes
tidyg2 <- tidyg %>% 
  activate(nodes) %>% 
  mutate(Connections = t_number[name], 
         Type = ifelse(name %in% org, "Microorganism", "Gene"),
         Sth = c(""))

p <- ggraph(tidyg2, layout = "lgl") +
  scale_edge_colour_gradient2() +
  geom_edge_link(aes(colour = Correlation)) +
  geom_node_point(aes(filter = Type != "Microorganism" & Connections >= 2, size = Connections)) +
  # geom_node_point(aes(filter = Type == "Microorganism", size = Connections)) +
  geom_node_label(aes(filter = Type == "Microorganism", label = name), 
                  size = 4, repel = TRUE) +
  scale_size(range = c(0.5, 2)) +
  theme_blank()
ggsave(paste0("Figures/",today,"_connections.png"))
  