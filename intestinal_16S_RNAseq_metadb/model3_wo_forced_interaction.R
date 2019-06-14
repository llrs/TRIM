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
C0 <- matrix(0, nrow = 5, ncol = 5)
# Use this to explore all models
C3 <- symm(C0, out_models[out_models$AVE_inner == max(out_models$AVE_inner), 
                                    grep("var.*", colnames(out_models))])
designs <- weight_design(weights = 11, size = 5, 
                         diff0 = which(lower.tri(C3) & C3 != 0))

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
saveRDS(out2, "model3_optimization_wo_forced_interaction.RDS")
out2 <- readRDS("model3_optimization_wo_forced_interaction.RDS")

C3 <- symm(designs[[1]], out2[out2$AVE_inner == max(out2$AVE_inner), grep("var.*", colnames(out2))])
dimnames(C3) <- list(names(A), names(A))
model3b3 <- sgcca(Ab, c1 = shrinkage, scheme = "centroid", C = C3, ncomp = rep(2, length(A)))
model3b3 <- improve.sgcca(model3b3, names(A))
saveRDS(model3b3, "model3_wo_forced_interaction.RDS")

l <- looIndex(size(A))
loo_model <- loo_functions(A, shrinkage)
result.out <- lapply(l, loo_model, model = C3)
saveRDS(result.out, "loo-model3_wo_forced_interaction.RDS")

pdf(paste0("Figures/", today, "_RGCCA_plots_model2.2.pdf"))

# Bootstrap of sgcca
boot <- boot_sgcca(A, C3, shrinkage, 1000)

saveRDS(boot, file = "bootstrap_model2.2.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

dev.off()