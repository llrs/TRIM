library("ggforce")
library("RGCCA")
library("Matrix")

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


points <- lapply(out, fuction(x) {
  for (o in x$Y) {
    for (i in x$X) {
      c(o, i)
    }
  }
})