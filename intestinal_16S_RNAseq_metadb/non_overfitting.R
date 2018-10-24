library("ggforce")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")


A <- readRDS("model3_TRIM.RDS")

# We cannnot comput eht tau.estimate for A[[1]]
shrinkage <- c(0.25670333, 0, 1, 1, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[[2]] <- tau.estimate(A[[2]])

# Experiment design for the complicated cases with too many computations possible to perform
designs <- weight_design(weights = 3, size = 5)
keep <- check_design(designs)

# Random subsample of 10% of the tryals
# Store all AVEs in the path so that it can be confirmed that it is the max value
s <- sample(designs[keep], size = min(length(designs[keep])*.1, 5000))
# Perform the sgcca on these samples
testing <- function(x, type, ...) {
  result.sgcca <- RGCCA::sgcca(C = x, 
                               scheme = type, 
                               verbose = FALSE, 
                               ...)
  analyze(result.sgcca)
}
# Estimated time of 8 hours
out <- sapply(s, testing, type = "centroid", A = A, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- t(out)
saveRDS(as.data.frame(out2), "sample_model3_boot.RDS")
# Linear model to see the weights and from there aim to the perfect AVE score
lmM0 <- lm(var12*var13*var14*var15*var23*var24*var25*var34*var35*var45~0+AVE_inner, data = as.data.frame(out2))
# Overfittin
lmM1 <- lm(AVE_inner~0+var12*var13*var14*var15*var23*var24*var25*var34*var35*var45, data = as.data.frame(out2))
library("broom")
library("dplyr")
tidy(lmM1) %>% 
  arrange(desc(abs(estimate)))
tidy(lmM0)
glance(lmM0)
glance(lmM1)

# Based on the top 4 of the random sample ie AVE_inner > 0.6
keep2 <- vapply(designs[keep], function(x) {
  x[4, 5] == 0 & x[1, 2] == 0 & x[2, 5] != 0 & x[2, 3] == 1 & x[3, 5] != 0
}, logical(1L))

out3 <- sapply(designs[keep][keep2], testing, 
               type = "centroid",
               A = A, c1 = shrinkage, 
               USE.NAMES = FALSE)
out3 <- t(out3)
saveRDS(as.data.frame(out3), "subset2_model3_boot.RDS")

# Select another range of 1000 of models and check which models are the best ones.
out <- predict(lmM0, data.frame(AVE_inner = seq(0.5, 1, by = 0.5)))
# Perhaps by looking at the best models and then for surrounding results with more affine weights
filter <- vapply(designs, function(x){
  l <- x[upper.tri(x)]
  all(l[c(3, 5, 8, 7, 9)] != 0) & all(l[c(1,2,4,6,10)] != 0)
}, FUN.VALUE = logical(1L))

out2 <- readRDS("sample_model3_boot.RDS")
columns <- grep("var", colnames(out3))
model3_best <- designs[[1]]
model3_best[upper.tri(model3_best)] <- unlist(out3[which.max(out3$AVE_inner), columns])
model3_best <- as.matrix(Matrix::forceSymmetric(model3_best))

# Leave one out procedure ####
# To asses if the selected model how well generalizes
l <- looIndex(size(A))
result.out <- lapply(l, function(x){
  
  RGCCA::sgcca(A = subsetData(A, x),
               C = model3_best, 
               scheme = "centroid", 
               verbose = FALSE, c1 = shrinkage
  )}) # SGCCA of the selected model leaving one sample each time out of order.
saveRDS(result.out, "loo-model3_best.RDS")
best <- sgcca(A, C = model3_best, c1 = shrinkage, verbose = FALSE, ncomp = rep(2, length(A)))
saveRDS(best, "model3_best.RDS")
# The results should be summarized/compared
# Which microorganisms are the same? Which genes are the same?
# Which AVE compared to original how well?
