library("ggforce")
library("RGCCA")
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

# Prepare the data:
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))

# Step 1 ####
# See the one with more AVE_inner on the boot with all models (but with just 3 weights)
out_models <- readRDS("sample_model3_boot.RDS")
fi <- out_models[out2$var12 != 0, ]

# Use this to explore all models
C3 <- symm(designs[[1]], fi[fi$AVE_inner == max(fi$AVE_inner), grep("var.*", colnames(fi))])
designs <- weight_design(weights = 11, size = 5, diff0 = which(lower.tri(C3) & C3 != 0))

# Perform the sgcca on these samples
testing <- function(x, type, ...) {
  result.sgcca <- RGCCA::sgcca(C = x, 
                                scheme = type, 
                                verbose = FALSE, 
                                scale = FALSE,
                                ...)
  analyze(result.sgcca)
}


# Explore for those cases
out <- sapply(designs, testing, type = "centroid", A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "sample_model3_forced_interaction.RDS")

# Select the best one
C3 <- symm(designs[[1]], out2[out2$AVE_inner == max(out2$AVE_inner), 
                              grep("var.*", colnames(out2))])
model3_best2 <- sgcca(Ab, c1 = shrinkage, scheme = "centroid", 
                      C = C3, ncomp = rep(2, length(Ab)))
model3_best2 <- improve.sgcca(model3_best2, names(Ab))
colnames(model3_best2$C) <- names(Ab)
rownames(model3_best2$C) <- names(Ab)
saveRDS(model3_best2, "model3_forced_interaction.RDS")

df <- data.frame(GE = model3_best2$Y[[1]][, 1],
                 M = model3_best2$Y[[2]][, 1],
                 D = model3_best2$Y[[3]][, 1],
                 L = model3_best2$Y[[4]][, 1],
                 Ti = model3_best2$Y[[5]][, 1])
df <- cbind(df, meta_r)
ggplot(df) +
  geom_point(aes(GE, M, color = Exact_location))
ggplot(df) +
  geom_point(aes(GE, M, color = IBD))
ggplot(df) +
  geom_point(aes(GE, M, color = SESCD_local))
