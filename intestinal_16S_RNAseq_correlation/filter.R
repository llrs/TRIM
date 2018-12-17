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

# Libraries ####
library("integration")
library("org.Hs.eg.db")
library("reactome.db")

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

select_genes_int <- function(file, expr) {
  # Load data from all the patients
  pre <- "../intestinal_16S_RNAseq_metadb"
  sgcca.centroid <- readRDS(file.path(pre, file))
  
  # Find outliers/important genes
  comp1 <- sgcca.centroid$a[[1]][, 1]
  outliers <- comp1 != 0
  comp1 <- comp1[outliers]
  
  keepGenes <- rownames(expr) %in% names(comp1)
  expr[keepGenes, ]
}

filter_sexual <- function(expr) {
  source("../genes_XY.R")
  sexual_related <- gsub("(.+)\\..*", "\\1", rownames(expr)) %in% 
  c(#bmX$ensembl_gene_id, 
    bmY$ensembl_gene_id)
  expr[!sexual_related, ]
}


readSGCCA <- function(file) {
  # Load data from all the patients
  sgcca.centroid <- readRDS(file)
  
  # Find outliers/important genes
  comp1 <- sgcca.centroid$a[[1]][, 1]
  outliers <- comp1 != 0
  comp1[outliers]
}

# Functions ####
filter_values <- function(comp1, cors, pval, threshold) {
  
  keepGenes <- colnames(cors) %in% names(comp1)
  cors <- cors[, keepGenes]
  pval <- pval[, keepGenes]
  
  cors <- cors[, !is.na(colnames(cors))]
  pval <- pval[, !is.na(colnames(pval))]
  
  message("Dimensions ", paste0(dim(cors), collapse = ", "))
  
  # Genes below the threshold
  keepCols <- apply(pval, 2, function(x){any(x < threshold)})
  keepRows <- apply(pval, 1, function(x){any(x < threshold)})
  
  keepCols[is.na(keepCols)] <- FALSE
  keepRows[is.na(keepRows)] <- FALSE
  
  if (sum(keepCols) == 0 || sum(keepRows) == 0) {
    stop("No relevant correlations with this threshold")
  }
  
  cors <- cors[keepRows, keepCols]
  pval <- pval[keepRows, keepCols]
  if (is.null(pval) || is.null(cors)) {
    stop("No relevant correlations with this threshold")
  }
  message("Dimensions ", paste0(dim(cors), collapse = ", "))
  
  list(cors = cors, pval = pval)
}

relevant <- function(comp1, cors, pval, threshold = 0.05) {
  l <- filter_values(comp1, cors, pval, threshold)
  pval <- l$pval
  cors <- l$cors
  if (ncol(pval) == 0) {
    stop("No relevant correlations with this threshold")
  }
  ind <- as.data.frame(which(pval < threshold, arr.ind = TRUE), 
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind)) # TODO test
  cor_pval <- apply(ind, 1, function(x){
    c("cors" = cors[x[1], x[2]],
      "pvalue" = pval[x[1], x[2]])
  })
  ind$row <- rownames(cors)[ind$row]
  ind$col <- colnames(cors)[ind$col]
  ind <- cbind(ind, t(cor_pval))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")
  
  ind <- ind[!duplicated(ind), ]
  ind <- ind[order(ind$Microorganism, ind$pvalue, decreasing = c(TRUE, FALSE)), ]
  rownames(ind) <- seq_len(nrow(ind))
  ind
}


sign_cor <- function(cors, pval, threshold = 0.05) {
  if (sum(pval < threshold, na.rm = TRUE) == 0) {
    stop("No relevant correlations with this threshold")
  }
  
  ind <- as.data.frame(which(pval < threshold, arr.ind = TRUE), 
                       stringAsFactors = FALSE)
  rownames(ind) <- seq_len(nrow(ind)) # TODO test
  cor_pval <- apply(ind, 1, function(x){
    c("cors" = cors[x[1], x[2]],
      "pvalue" = pval[x[1], x[2]])
  })
  ind$row <- rownames(cors)[ind$row]
  ind$col <- colnames(cors)[ind$col]
  ind <- cbind(ind, t(cor_pval))
  colnames(ind) <- c("Microorganism", "Gene", "Correlation", "pvalue")
  
  ind <- ind[!duplicated(ind), ]
  ind <- ind[order(ind$Microorganism, ind$pvalue, decreasing = c(TRUE, FALSE)), ]
  rownames(ind) <- seq_len(nrow(ind))
  ind
}
# Expects genes in rows and species at the columns
plot_cor <- function(file, cors, pval, threshold, label) {
  l <- filter_values(file, cors, pval, threshold)
  cors <- l$cors
  
  cors <- cors[!duplicated(rownames(cors)), ]
  colors_g <- ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                            midpoint = 0, limits = c(-1, 1))
  heatmaply(cors, name = "Cor",
            ylab = "Genes",
            xlab = "Genus",
            scale_fill_gradient_fun = colors_g,
            file = paste0("Figures/", today, "heatmap", label,".html"))
  
}


ensembl2symbol <- function(x) {
  all_samples_symbol <- x
  all_samples_symbol$Gene <- trimVer(all_samples_symbol$Gene)
  all_samples_symbol$Gene <- mapIds(org.Hs.eg.db, keys = all_samples_symbol$Gene,
                                    keytype = "ENSEMBL", column = "SYMBOL")
  all_samples_symbol <- all_samples_symbol[!is.na(all_samples_symbol$Gene), ]
  all_samples_symbol[!duplicated(all_samples_symbol), ]
}
write_cor <- function(x, file){
  all_samples_symbol <- ensembl2symbol(x)
  write.csv(all_samples_symbol, file = file, row.names = FALSE, na = "")
}

#' Enrichment by microorganisms
#' 
#' Function to split by microorganism and compare if they are enriched in some 
#' pathways
#' @param all All the genes that are found relevant (weight != 0)
#' @param x The table of microoganisms, genes, correlations and adjusted p-values.
#' @return A table with the pathways enriched for every microorganism
pathsPerMicro <- function(x, all_genes){
  genesOrig <- trimVer(names(all_genes))
  x[, "Gene"] <- trimVer(x[, "Gene"])
  perMicro <- split(x, x[, "Microorganism"])
  entrezID <- mapIds(org.Hs.eg.db, keys = genesOrig, keytype = "ENSEMBL", 
                     column = "ENTREZID")
  entrezID <- entrezID[!is.na(entrezID)]
  paths2genes <- access_reactome()
  genes <- unlist(paths2genes, use.names = FALSE)
  pathways <- rep(names(paths2genes), lengths(paths2genes))
  message("Calculating the enrichment")
  T2G <- cbind.data.frame(pathways, genes)
  lapply(perMicro, function(x) {
    significant <- x$Gene
    relevant <- entrezID[significant]
    relevant <- relevant[!is.na(relevant)]
    if (length(relevant) <= 10) {
      return(NULL)
    }
    enrich <- clusterProfiler::enricher(gene = relevant , 
                                        # universe = entrezID, 
                                        TERM2GENE = T2G)
    as.data.frame(enrich)
    enrich <- as.data.frame(enrich)
    if (nrow(enrich) >= 1) {
      enrich$Description <- mapIds(reactome.db, keys = rownames(enrich), 
                                   keytype = "PATHID", column = "PATHNAME")
    }
    enrich
  })
}

# Output ####
library("reactome.db")
# #  The numbers come from the number of samples in each correlation
threshold <- 0.05
all_comp <- readSGCCA("../intestinal_16S_RNAseq_integration/sgcca.RDS")
all_samples_ensembl <- relevant(all_comp, all_s, pall_s, threshold)
o2 <- pathsPerMicro(all_samples_ensembl, all_comp)
all_samples_ensembl_cor <- sign_cor(all_s, pall_s, threshold)
write_cor(all_samples_ensembl, file = paste0(today, "_", type, "_correlation_all_model0.csv"))
# write_cor(all_samples_ensembl_cor, file = paste0(today, "_sign_", type, "_correlation_all.csv"))
lapply(names(o2), function(x) {
  if (nrow(o2[[x]]) > 1 && !is.null(o2[[x]])) {
    write.csv(o2[[x]], row.names = FALSE, na = "", 
              file = paste0(x, "_pathways_by_cor_all.csv"))
  }
})
colon_samples_ensembl <- relevant(all_comp, colon, pcolon, threshold)
o2 <- pathsPerMicro(colon_samples_ensembl, all_comp)
colon_samples_ensembl_cor <- sign_cor(colon, pcolon, threshold)
write_cor(colon_samples_ensembl, file = paste0(today, "_", type, "_correlation_colon_model0.csv"))

ileum_samples_ensembl <- relevant(all_comp, ileum, pileum, threshold)
o2 <- pathsPerMicro(ileum_samples_ensembl, all_comp)
lapply(names(o2), function(x) {
  if (nrow(o2[[x]]) > 1 && !is.null(o2[[x]])) {
    write.csv(o2[[x]], row.names = FALSE, na = "", 
              file = paste0(x, "_pathways_by_cor_ileum.csv"))
  }
})
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

# Bootstrapping to asses how probable is to have such enrichment in those p-values matrix
boots_corr <- function(pvalues, component, iter = 10000) {
  b <- vapply(seq_len(iter), function(x){
    sel <- sample(x = seq_len(ncol(pvalues)), size = length(component))
    sum(pvalues[, sel] < 0.05, na.rm = TRUE)
  }, numeric(1L))
  n <- sum(pvalues[, colnames(pvalues) %in% names(component)] < 0.05, na.rm = TRUE)
  sum(b >= n)/length(b)
}

c_disease <- boots_corr(pdisease, disease_comp)
c_all <- boots_corr(pall_s, all_comp)
c_controls <- boots_corr(pcontrols, contr_comp)


plot_single_cor <- function(x, gene, y, microorganism, colr, case, cor_val) {
  genes <- trimVer(gene)
  symbol <- tryCatch({mapIds(
    org.Hs.eg.db, key = genes, keytype = "ENSEMBL",
    column = "SYMBOL"
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

library("ggplot2")
library("ggraph")
library("tidygraph")
library("ggnetwork")
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
  