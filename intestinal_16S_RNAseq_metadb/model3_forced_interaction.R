library("ggforce")
library("RGCCA2")
library("integration")
library("fgsea")
library("dplyr")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")

A <- readRDS("model3_TRIM.RDS")
meta_r <- readRDS( "meta.RDS")

# We cannnot comput eht tau.estimate for A[[1]] calculated on server
# original value 0.25670333
# last calculated value is of 0.249488046688595
shrinkage <- c(0.249488046688595, 0, 1, 1, 1) 
shrinkage[[2]] <- tau.estimate(A[[2]]) # 0.286506412433534

# Experiment design for the complicated cases 
# with too many computations possible to perform it only tests for 3 weighs per edge
designs <- weight_design(weights = 11, size = 5, c(2, 4, 8, 9, 15))
keep <- vapply(designs, RGCCA2::correct, logical(1L))
designs <- designs[keep]

keep2 <- vapply(designs, function(x){
  x[2, 1] != 0 && x[2, 4] != 0 && x[1, 4] != 0 && x[3, 5] != 0 && x[2, 3] != 0
}, logical(1L))

designs <- designs[keep2]

# Step 1 ####
# Random subsample of 10% of the trials
# Store all AVEs in the path so that it can be confirmed that it is the max value
set.seed(4672679)
s <- sample(designs, size = min(length(designs)*.1, 10000))
# Perform the sgcca on these samples
testing <- function(x, type, ...) {
  result.sgcca <- RGCCA2::sgcca(C = x, 
                                scheme = type, 
                                verbose = FALSE, 
                                scale = FALSE,
                                ...)
  analyze(result.sgcca)
}
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# Estimated time of 8 hours
out <- sapply(s, testing, type = "centroid", A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "sample_model3_forced_interaction.RDS")

out2 %>% 
  top_n(5, AVE_inner) %>% 
  select(AVE_inner, AVE_outer, var12, var13, var23, 
         var14, var24, var34, var15, var25, var35, var45) %>% 
  arrange(desc(AVE_inner))
stop("Visual inspection of the top 5")
# step 2 ####

C3 <- symm(designs[[1]], out2[out2$AVE_inner == max(out2$AVE_inner), grep("var.*", colnames(out2))])
model3_best2 <- sgcca(Ab, c1 = shrinkage, scheme = "centroid", C = C3)
saveRDS(model3_best2, "model3_forced_interaction.RDS")

C3 <- subSymm(C3, 1, 2, 0)
model3b3 <- sgcca(Ab, c1 = shrinkage, scheme = "centroid", C = C3)
saveRDS(model3b3, "model3_wo_forced_interaction.RDS")

df <- data.frame(GE = model3b3$Y[[1]][, 1],
                 M = model3b3$Y[[2]][, 1],
                 D = model3b3$Y[[3]][, 1],
                 L = model3b3$Y[[4]][, 1],
                 Ti = model3b3$Y[[5]][, 1])
df <- cbind(df, meta_r)
ggplot(df) +
  geom_point(aes(GE, M, color = Exact_location))
ggplot(df) +
  geom_point(aes(GE, M, color = IBD))
ggplot(df) +
  geom_point(aes(GE, M, color = SESCD_local))
