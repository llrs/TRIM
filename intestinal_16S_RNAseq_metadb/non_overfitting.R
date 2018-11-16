library("ggforce")
library("RGCCA")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")
library("fgsea")


A <- readRDS("model3_TRIM.RDS")

# We cannnot comput eht tau.estimate for A[[1]] calculated on server
shrinkage <- c(0.25670333, 0, 1, 1, 1) 
shrinkage[[2]] <- tau.estimate(A[[2]])

# Experiment design for the complicated cases 
# with too many computations possible to perform it only tests for 3 weighs per edge
designs <- weight_design(weights = 3, size = 5)
keep <- vapply(designs, correct, logical(1L))
designs <- designs[keep]

# Step 1 ####
# Random subsample of 10% of the trials
# Store all AVEs in the path so that it can be confirmed that it is the max value
set.seed(4672679)
s <- sample(designs, size = min(length(designs)*.1, 1000))
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
out2 <- as.data.frame(t(out))
saveRDS(out2, "sample_model3_boot.RDS")
library("dplyr")
out2 %>% 
  top_n(5, AVE_inner) %>% 
  select(AVE_inner, AVE_outer, var12, var13, var23, 
         var14, var24, var34, var15, var25, var35, var45) %>% 
  arrange(desc(AVE_inner))
stop("Visual inspection of the top 5")
# step 2 ####

# Based on the top 5 of the random sample ie AVE_inner > 0.6
# Select another range of 37000 of models and check which models are the best ones.
best_keep <- vapply(designs, function(x) {
  x[4, 5] == 0 & x[3, 4] == 0 & x[3, 5] != 0
}, logical(1L))

out3 <- sapply(designs[best_keep], testing, 
               type = "centroid",
               A = A, c1 = shrinkage, 
               USE.NAMES = FALSE)
out3 <- as.data.frame(t(out3))
saveRDS(out3, "subset2_model3_boot.RDS")
stop("Think again!")
# Step 3 ####
out2 <- readRDS("sample_model3_boot.RDS")
out3 <- readRDS("subset2_model3_boot.RDS")
out <- rbind(out2, out3)
columns <- grep("var", colnames(out))
out <- out[!duplicated(out), ]
model3_best <- designs[[1]]
model3_best <- symm(model3_best, out[which.max(out$AVE_inner), columns])
colnames(model3_best) <- names(A)
rownames(model3_best) <- names(A)
out4 <- out
out4$weights <- as.factor(out4$weights)
ggplot(out4) + geom_boxplot(aes(weights, AVE_inner))
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

model3_best2 <- subSymm(model3_best, 1, 1, 1)
best_interaction <- sgcca(A, C = model3_best2, c1 = shrinkage, verbose = FALSE, ncomp = rep(2, length(A)))
saveRDS(best_interaction, "model3_best_interaction.RDS")

result.i <- lapply(l, function(x){
  
  RGCCA::sgcca(A = subsetData(A, x),
               C = model3_best2, 
               scheme = "centroid", 
               verbose = FALSE, c1 = shrinkage
  )}) # SGCCA of the selected model leaving one sample each time out of order.
saveRDS(result.i, "loo-model3_best_interaction.RDS")
a <- best$a[[1]][, 1]
ai <-  best_interaction$a[[1]][, 1]
hist(a[a != 0 | ai != 0] - ai[a != 0 | ai != 0]) # Normal distribution: randomness