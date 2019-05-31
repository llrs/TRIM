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

# Experiment design for the complicated cases 
# with too many computations possible to perform it only tests for 3 weighs per edge
designs <- weight_design(weights = 11, size = 5, c(4, 8, 9, 15))
keep <- vapply(designs, RGCCA::correct, logical(1L))
designs <- designs[keep]

keep2 <- vapply(designs, function(x){
  x[2, 1] == 0 && x[2, 4] != 0 && x[1, 4] != 0 && x[3, 5] != 0 && x[2, 3] != 0
}, logical(1L))

designs <- designs[keep2]

# Step 1 ####
# Perform the sgcca on these samples
testing <- function(x, type, ...) {
  tryCatch({
  result.sgcca <- RGCCA::sgcca(C = x, 
                                scheme = type, 
                                verbose = FALSE, 
                                scale = FALSE,
                                ...)
  analyze(result.sgcca)}, error = function(e){
    out <- rep(NA, 24)
    names(out) <- c(
      "vs12", "vs13", "vs23", "vs14", "vs24", "vs34", "vs15",
      "vs25", "vs35", "vs45", "AVE_inner", "AVE_outer", "cc1", "var12", "var13",
      "var23", "var14", "var24", "var34", "var15", "var25", "var35",
      "var45", "weights"
    )
  }
    
  )
}
Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# Estimated time of 8 hours
out <- sapply(designs, testing, type = "centroid", A = Ab, c1 = shrinkage, USE.NAMES = FALSE)
out2 <- as.data.frame(t(out))
saveRDS(out2, "model3_optimization_wo_forced_interaction.RDS")
out2 <- readRDS("model3_optimization_wo_forced_interaction.RDS")
# This code is needed to remove the factors introduced by handling the error :\
keep <- out2$vs12 != "vs12"
out2 <- out2[keep, ]
out2 <- sapply(out2, function(x){
  as.numeric(levels(x))[x]
})
out2 <- as.data.frame(out2)

out2 %>% 
  select(AVE_inner, AVE_outer, matches("var")) %>% 
  arrange(desc(AVE_inner)) %>% 
  head()
stop("Visual inspection of the top 5")

C3 <- symm(designs[[1]], out2[out2$AVE_inner == max(out2$AVE_inner), grep("var.*", colnames(out2))])
dimnames(C3) <- list(names(A), names(A))
model3b3 <- sgcca(Ab, c1 = shrinkage, scheme = "centroid", C = C3)
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