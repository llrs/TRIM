library("ggforce")
library("RGCCA")
library("BiocParallel")
library("integration")
library("fgsea")
library("broom")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

# Load data
otus_table_i <- readRDS("../intestinal_16S_RNAseq_integration/otus_table_norm_RNAseq.RDS")
otus_tax_i <- readRDS("../intestinal_16S_RNAseq_metadb/otus_tax.RDS")
expr <- readRDS("../intestinal_16S_RNAseq_metadb/expr.RDS")
meta_r <- readRDS("../intestinal_16S_RNAseq_metadb/meta.RDS")

keep <- meta_r$Exact_location != "ILEUM" | is.na(meta_r$Exact_location)

otus_table_i <- otus_table_i[, keep]
expr <- expr[, keep]
metadb <- meta_r[keep, ]

A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), meta = metadb)
A[1:2] <- clean_unvariable(A[1:2])
saveRDS(A, "colon.RDS")
# shrinkage <- sapply(A, tau.estimate)
shrinkage <- c(RNAseq = 0.343325970491348, `16S` = 0.424174371340925)
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0 <- subSymm(C[1:2, 1:2], "16S", "RNAseq", 1)
m0 <- sgcca(A[1:2], model0, c1 = shrinkage, ncomp = rep(2, 2))
m0 <- improve.sgcca(m0, names(A)[1:2])
saveRDS(m0, "model0.RDS")

m0$AVE
metadata <- model_RGCCA(metadb, c("AgeDiag", "AGE_SAMPLE", "Transplant",
                                  "ID","SEX", "Surgery", "Treatment"))
metadata$AgeDiag[is.na(metadata$AgeDiag)] <- 0 # Time has NA values

model1 <- subSymm(C, "16S", "meta", 1)
model1 <- subSymm(model1, "meta", "RNAseq", 1)
A$meta <- metadata
shrinkage1 <- c(shrinkage, 1)
m1 <- sgcca(A, model1, c1 = shrinkage1, ncomp = rep(2, 3))
m1 <- improve.sgcca(m1, names(A))
saveRDS(m1, "model1.RDS")

model1.1 <- subSymm(model1, "16S", "RNAseq", 1)
m1.1 <- sgcca(A, model1.1, c1 = shrinkage1, ncomp = rep(2, 3))
m1.1 <- improve.sgcca(m1.1, names(A))
saveRDS(m1.1, "model1.1.RDS")

designs1.2 <- weight_design(weights = 11, size = 3)
k <- vapply(designs1.2, correct, logical(1L))
designs1.2 <- designs1.2[k]

Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

testing <- function(x, ...) {
  result.sgcca <- RGCCA::sgcca(C = x, 
                               verbose = FALSE, 
                               scale = FALSE,
                               ...)
  analyze(result.sgcca)
}

out <- sapply(designs1.2, testing, A = Ab, c1 = shrinkage1, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model1.2_testing.RDS")


columns <- grep("var", colnames(out))
model1.2 <- symm(model1.1, out[which.max(out$AVE_inner), columns])
m1.2 <- sgcca(A, model1.2, c1 = shrinkage1, ncomp = rep(2, 3))
m1.2 <- improve.sgcca(m1.2, names(A))
saveRDS(m1.2, "model1.2.RDS")

Time <- model_RGCCA(metadb, c("AgeDiag", "AGE_SAMPLE", "Transplant"))
Demographics <- model_RGCCA(metadb, c("ID","SEX", "Surgery", "Treatment"))
Time$AgeDiag[is.na(Time$AgeDiag)] <- 0 # Time has NA values

A <- c(A[1:2], list(Time = Time, Demographics = Demographics))
shrinkage2 <- c(shrinkage, 1, 1)
names(shrinkage2) <- names(A)
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model2 <- subSymm(C, "16S", "RNAseq", 1)
model2 <- subSymm(model2, "RNAseq", "Demographics", 1)
model2 <- subSymm(model2, "RNAseq", "Time", 1)
model2 <- subSymm(model2, "16S", "Demographics", 1)
model2 <- subSymm(model2, "16S", "Time", 1)
m2 <- sgcca(A, model2, c1 = shrinkage2, ncomp = rep(2, length(A)))
m2 <- improve.sgcca(m2, names(A))
saveRDS(m2, "model2.RDS")

Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
designs2.1 <- weight_design(weights = 3, size = 4)
k <- vapply(designs2.1, correct, logical(1L))
designs2.1 <- designs2.1[k]
out <- sapply(designs2.1, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model2.1_testing.RDS")

columns <- grep("var", colnames(out))
model2.1 <- symm(model2, out[which.max(out$AVE_inner), columns])
m2.1 <- sgcca(A, model2.1, c1 = shrinkage2, ncomp = rep(2, 4))
m2.1 <- improve.sgcca(m2.1, names(A))
saveRDS(m2.1, "model2.1.RDS")

designs2.2 <- weight_design(weights = 11, size = 4, 
                            diff0 = which(lower.tri(model2.1) & model2.1 != 0))
out <- sapply(designs2.2, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model2.2_testing.RDS")

model2.2 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
m2.2 <- sgcca(A, model2.2, c1 = shrinkage2, ncomp = rep(2, 4))
m2.2 <- improve.sgcca(m2.2, names(A))
saveRDS(m2.2, "model2.2.RDS")


columns <- grep("var", colnames(out))
model2.3 <- symm(model2, out[out$var12 != 0, ][which.max(out$AVE_inner), columns])
designs2.3 <- weight_design(weights = 11, size = 4, 
                            diff0 = which(lower.tri(model2.3) & model2.3 != 0))
out <- sapply(designs2.3, testing, A = Ab, c1 = shrinkage2, USE.NAMES = FALSE)
out <- as.data.frame(t(out))
saveRDS(out, "model2.3_testing.RDS")


model2.3 <- symm(C, out[which.max(out$AVE_inner), grep("var", colnames(out))])
m2.3 <- sgcca(A, model2.3, c1 = shrinkage2, ncomp = rep(2, 4))
m2.3 <- improve.sgcca(m2.3, names(A))
saveRDS(m2.3, "model2.3.RDS")


a <- mget(ls(pattern = "^m[0-9]"))
sapply(a, function(x){x$AVE$AVE_inner})