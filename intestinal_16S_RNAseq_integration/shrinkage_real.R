library("RGCCA")
library("integration")
library("ggplot2")
library("patchwork")

# Save
otus_table_i <- readRDS("otus_table.RDS")
expr <- readRDS("expr.RDS")
meta_r <- readRDS( "meta.RDS")

# Read
A <- readRDS(file = "TRIM.RDS")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model0 <- subSymm(C, "16S", "RNAseq", 1)

# We cannnot comput the tau.estimate for A[[1]]
min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
})

shrinkage <- expand.grid(
  seq(from = min_shrinkage[1], to = 1, length.out = 10),
  seq(from = min_shrinkage[2], to = 1, length.out = 10))

# The one tau.estimate says:
shrinkage_opt <- c(0.249488046688595, 0) 
shrinkage_opt[2] <- tau.estimate(A[[2]])
shrinkage_opt2 <- shrinkage_opt
shrinkage_opt2[1] <- 0.1

shrinkage <- rbind(shrinkage, shrinkage_opt, shrinkage_opt2)

testing <- function(x, type, ...) {
  result.sgcca <- RGCCA::sgcca(c1 = x, 
                               scheme = type, 
                               verbose = FALSE, 
                               scale = FALSE,
                               ...)
  selected <- vapply(result.sgcca$a, function(y){
    sum(y[, 1] != 0)
  }, numeric(1L))
  
  c("AVE_inner" = result.sgcca$AVE$AVE_inner, 
    "AVE_outer" = result.sgcca$AVE$AVE_outer, 
    selected)
}

Ab <- lapply(A, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
# Estimated time of 1 hours
bench::system_time(out <- apply(shrinkage, 1, testing, 
                         type = "centroid", A = Ab, C = model0))
o <- t(out)
beepr::beep()

dat <- cbind.data.frame(o, shrinkage)
saveRDS(dat, "shrinkage_real.RDS")
r <- ggplot(dat) +
  geom_point(aes(Var1, RNAseq, color = Var2)) +
  scale_color_viridis_c() +
  labs(x = "Shrinkage RNAseq", color = "Shrinkage 16S", title = "Effect of shrinkage parameter", 
       y = "# genes") +
  geom_vline(xintercept = 0.1)
s <- ggplot(dat) +
  geom_point(aes(Var2, `16S`, color = Var1)) +
  scale_color_viridis_c() +
  labs(color = "Shrinkage RNAseq", x = "Shrinkage 16S", title = "Effect of shrinkage parameter", 
       y = "# OTUs")
r + s
r2 <- ggplot(dat) +
  geom_point(aes(RNAseq, AVE_inner,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_inner", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# genes", size = "Srhinkage 16S") +
  scale_size_continuous()
r3 <- ggplot(dat) +
  geom_point(aes(`16S`, AVE_inner,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_inner", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# 16S", size = "Srhinkage 16S") +
  scale_size_continuous()
s2 <- ggplot(dat) +
  geom_point(aes(RNAseq, AVE_outer,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_outer", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# genes", size = "Srhinkage 16S") +
  scale_size_continuous()
r2 + s2

ggplot(dat) + 
  geom_point(aes(Var1, AVE_inner, color = Var2, size = log10(RNAseq)))