library("integration")
library("dplyr")
library("BaseSet")
library("ggplot2")
library("org.Hs.eg.db")
library("clusterProfiler")  
library("ReactomePA")

# Taxonomy
tax <- readRDS("intestinal_16S_RNAseq_metadb/otus_tax.RDS")

# Evaluate all models ###
# model 0
folder0 <- "intestinal_16S_RNAseq_integration"
model0 <- readRDS(file.path(folder0, "sgcca.RDS"))
model0_loo <- readRDS(file.path(folder0, "loo-model0.RDS"))

# model 1 without interaction
folder1 <- "intestinal_16S_RNAseq_metadb"
# model 2 without interaction
model1.2 <- readRDS(file.path(folder1, "model2_best.RDS"))
model1.2_loo <- readRDS(file.path(folder1, "loo-model2_best.RDS"))
# model 3
model2.2 <- readRDS(file.path(folder1, "model3_wo_forced_interaction.RDS"))
model2.2_loo <- readRDS(file.path(folder1, "loo-model3_wo_forced_interaction.RDS"))

relations_sgcca2 <- function(x, nm) {
  elements <- lapply(x$a, function(x){rownames(x)[x[, 1] != 0]})
  elements <- unlist(elements[1:2], use.names = FALSE)
  sets <- rep(nm, length(elements))
  data.frame(elements = elements, sets = sets)
}

fuzzy_model_set <- function(model_orig, model_loo, name) {
  l <- vector("list", length = length(model_loo))
  
  for (i in seq_along(model_loo)) {
    l[[i]] <- relations_sgcca2(model_loo[[i]], name)
  }
  
  orig <- relations_sgcca2(model_orig, name)
  r <- do.call(rbind, l)
  r <- rbind(r, orig)
  
  r %>% 
    group_by(elements, sets) %>% 
    count() %>% 
    mutate(fuzzy = n/(length(model_loo) + 1 )) %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    dplyr::select(elements, sets, fuzzy)
}

# Create the data.frame with the relations probability and the strength on the original data
f0 <- fuzzy_model_set(model0, model0_loo, "0")
a <- rbind(model0$a[[1]][, 1, drop = FALSE], model0$a[[2]][, 1, drop = FALSE])
f0 <- merge(f0, a, by.x = "elements", by.y = "row.names", all.x = TRUE, all.y = FALSE)
f1.2 <- fuzzy_model_set(model1.2, model1.2_loo, "1.2")
a <- rbind(model1.2$a[[1]][, 1, drop = FALSE], model1.2$a[[2]][, 1, drop = FALSE])
f1.2 <- merge(f1.2, a, by.x = "elements", by.y = "row.names", all.x = TRUE, all.y = FALSE)
f2.2 <- fuzzy_model_set(model2.2, model2.2_loo, "2.2")
a <- rbind(model2.2$a[[1]][, 1, drop = FALSE], model2.2$a[[2]][, 1, drop = FALSE])
f2.2 <- merge(f2.2, a, by.x = "elements", by.y = "row.names", all.x = TRUE, all.y = FALSE)

relations <- rbind(f0, f1.2, f2.2)


models <- tidySet(relations) %>% 
  mutate_element(Type = if_else(grepl("^ENSG", elements), "Gene", "OTU"))

elements_model0 <- models %>% 
  filter_set(sets == "0") %>% 
  elements() %>% 
  pull(elements) %>% 
  as.character()

third <- intersection(models, c("2.2", "1.2"), keep = TRUE, name = "X.2") %>% 
  subtract(set_in = "X.2", not_in = "0", name = "CD", keep = TRUE) %>% 
  intersection(sets = c("0", "2.2", "1.2"), name = "Baseline", keep = TRUE)

third %>% 
  relations() %>% 
  group_by(sets) %>% 
  count()

# Which models have more elements that remain constant thorough the LOO ?
models %>% 
  relations() %>% 
  group_by(sets, fuzzy) %>% 
  count() %>% 
  group_by(sets) %>% 
  mutate(total = sum(n), freq = n / total) %>% 
  ggplot() + 
  geom_point(aes(fuzzy, freq, col = sets), alpha = 0.7) +
  ylim(c(0, 1)) +
  NULL

third %>% 
  filter_relation(fuzzy == 1) %>%
  relations() %>% 
  group_by(sets) %>% 
  count()

m <- models %>% 
  filter_relation(fuzzy == 1) %>%
  incidence()
m[m != 0] <- 1
m <- m[, c(2, 3, 1)]
UpSetR::upset(as.data.frame(m), order.by = "freq", keep.order = TRUE,
              sets = rev(c("0", "1.2" ,"2.2")), sets.x.label = "Model size")

# Genes function
third <- third %>%
  mutate_element(Ensembl = trimVer(elements),
                 Symbol = mapIds(org.Hs.eg.db, keys = Ensembl, 
                                 keytype = "ENSEMBL", 
                                 column = "SYMBOL", multiVals = "first"),
                 entrez = mapIds(org.Hs.eg.db, keys = Ensembl,
                                 keytype = "ENSEMBL", 
                                 column = "ENTREZID", multiVals = "first"))

o <- third %>% 
  filter(sets %in% c("CD", "Baseline")) %>% 
  as.data.frame()
  # dplyr::select(-sets)
  # filter(fuzzy > 0.95) %>% 
  NA
  
exclusive <- third %>% 
  subtract(set_in = "1.2", not_in = c("0", "2.2"), name = "exclusive 1.2", keep = TRUE) %>% 
  subtract(set_in = "2.2", not_in = c("0", "1.2"), name = "exclusive 2.2", keep = TRUE) %>% 
  subtract(set_in = "0", not_in = c("1.2", "2.2"), name = "exclusive 0", keep = TRUE) %>% 
  filter_set(sets %in% c("exclusive 1.2", "exclusive 2.2", "exclusive 0")) %>% 
  as.data.frame() %>% 
  arrange(desc(fuzzy), Symbol, desc(Type))

taxes <- as.data.frame(tax) %>% 
  tibble::rownames_to_column()
exclusive_full <- exclusive %>% 
  left_join(taxes, by = c("elements" = "rowname"))

write.csv(exclusive_full[exclusive_full$sets == "exclusive 1.2", ], na = "", 
          file = "exclusive_model_1.2.csv", row.names = FALSE)
write.csv(exclusive_full[exclusive_full$sets == "exclusive 2.2", ], na = "", 
          file = "exclusive_model_2.2.csv", row.names = FALSE)
subset_third <- third %>% 
  filter(sets %in% c("CD", "Baseline"))

write.csv(as.data.frame(o)[, -2], file = "relevant_genes_different_models.csv", 
          row.names = FALSE)

universe <- mapIds(org.Hs.eg.db, keys = trimVer(rownames(model0$a$RNAseq)), 
                   keytype = "ENSEMBL", column = "ENTREZID")

x <- enrichPathway(gene = exclusive_full[exclusive_full$sets == "exclusive 2.2", "entrez"], pvalueCutoff = 0.05, readable = TRUE,
                   universe = universe)
barplot(x, showCategory = 30)
dotplot(x, showCategory = 30)
emapplot(x, showCategory = 30)
cnetplot(x, categorySize = "pvalue")

# Taxonomic

OTUs <- third %>% 
  filter_relation(fuzzy == 1) %>% 
  filter_set(sets %in% c("CD", "Baseline")) %>% 
  elements() %>%
  filter(Type == "OTU") %>% 
  pull(elements) %>% 
  trimVer()
table(tax[OTUs, "Family"])

third %>% 
  filter_relation(fuzzy == 1) %>% 
  filter_set(sets %in% c("CD", "Baseline")) %>% 
  elements() %>%
  group_by(Type) %>% 
  count()
