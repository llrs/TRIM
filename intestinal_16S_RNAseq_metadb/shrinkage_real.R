library("RGCCA2")
library("integration")
library("ggplot2")
library("patchwork")
library("dplyr")

# Read data
A <- readRDS(file = "TRIM.RDS")

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
model1 <- subSymm(C, "16S", "metadata", 1)
model1 <- subSymm(model1, "RNAseq", "metadata", 1)
model2 <- subSymm(model1, "RNAseq", "16S", 1)
model2b <- subSymm(model1, "RNAseq", "metadata", 0.1)

A2 <- readRDS("model3_TRIM.RDS")
C2 <- matrix(
  0, ncol = length(A2), nrow = length(A2),
  dimnames = list(names(A2), names(A2))
)
model3 <- subSymm(C2, "16S", "RNAseq", 1)
model3 <- subSymm(model3, "16S", "Demographics", 1)
model3 <- subSymm(model3, "16S", "Location", 1)
model3 <- subSymm(model3, "16S", "Time", 1)
model3 <- subSymm(model3, "RNAseq", "Demographics", 1)
model3 <- subSymm(model3, "RNAseq", "Location", 1)
model3 <- subSymm(model3, "RNAseq", "Time", 1)

model3b <- subSymm(C2, "RNAseq", "Location", 1)
model3b <- subSymm(model3b, "16S", "Location", 0.5)
model3b <- subSymm(model3b, "16S", "Demographics", 1)
model3b <- subSymm(model3b, "Time", "Demographics", 1)

# We cannnot comput the tau.estimate for A[[1]]
min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
})

shrinkage <- expand.grid(
  seq(from = min_shrinkage[1], to = 1, length.out = 10),
  seq(from = min_shrinkage[2], to = 1, length.out = 10),
  1)

# The one tau.estimate says:
shrinkage_opt <- c(0.249488046688595, 0, 1) 
shrinkage_opt[2] <- tau.estimate(A[[2]])
shrinkage_opt2 <- shrinkage_opt
shrinkage_opt2[1] <- 0.1

shrinkage <- rbind(shrinkage, shrinkage_opt, shrinkage_opt2)

# It uses my own improved modification of RGCCA which is faster
# However it 
testing <- function(x, ...) {
  result.sgcca <- RGCCA2::sgcca(c1 = x,  
                               scheme = "centroid", 
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
out <- apply(shrinkage, 1, testing, A = Ab, C = model1)
o <- t(out)
dat <- cbind.data.frame(o, shrinkage)
saveRDS(dat, "shrinkage_model1.RDS")

out <- apply(shrinkage, 1, testing, A = Ab, C = model2)
o <- t(out)
dat <- cbind.data.frame(o, shrinkage)
saveRDS(dat, "shrinkage_model2.RDS")

out <- apply(shrinkage, 1, testing, A = Ab, C = model2b)
o <- t(out)
dat <- cbind.data.frame(o, shrinkage)
saveRDS(dat, "shrinkage_model2b.RDS")

Ab2 <- lapply(A2, function(x) scale2(x, bias = TRUE)/sqrt(NCOL(x)))
shrinkage2 <- cbind(shrinkage, "Var4" = 1, "Var5" = 1)
out <- apply(shrinkage2, 1, testing, A = Ab2, C = model3)
o <- t(out)
dat <- cbind.data.frame(o, shrinkage2)
saveRDS(dat, "shrinkage_model3.RDS")


out <- apply(shrinkage2, 1, testing, A = Ab2, C = model3b)
o <- t(out)
dat <- cbind.data.frame(o, shrinkage2)
saveRDS(dat, "shrinkage_model3b.RDS")

theme_set(theme_bw())
plot_shrinkage <- function(real, model) {
  
  q <- ggplot(real) + 
    geom_point(aes(Var1,  Var2, color = AVE_inner, alpha = AVE_inner)) +
    scale_color_viridis_c(begin = 0, end = 1) +
    ggtitle(model) +
    labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
         size = "# genes") +
    guides(alpha = FALSE) +
    scale_alpha(range = c(0, 1))
  print(q)
  p <- ggplot(real) + 
    geom_point(aes(Var1,  Var2, color = AVE_outer, alpha = AVE_inner)) +
    scale_color_viridis_c(begin = 0, end = 1) +
    ggtitle(model) +
    labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
         size = "# genes") +
    guides(alpha = FALSE) +
    scale_alpha(range = c(0, 1))
  c(p, q)
}

flies <- list.files(pattern = "shrinkage_model.*.RDS")
flies <- flies[c(1, 3, 2, 5, 4)]
helper <- function(x){
  model <- strsplit(x, "_|\\.")[[1]][2]
  out <- readRDS(x)
  cbind.data.frame(out[, c("AVE_inner", "AVE_outer", "Var1", "Var2")], 
                   model = model)
}
o <- lapply(flies, helper)
a <- Reduce(rbind.data.frame, o)

o0 <- readRDS("../intestinal_16S_RNAseq_integration/shrinkage_real.RDS")
a2 <- rbind.data.frame(cbind.data.frame(o0[, c("AVE_inner", "AVE_outer", "Var1", "Var2")], 
                          model = "0"), a)
levels(a2$model)[levels(a2$model) == "model1"] <- "1"
levels(a2$model)[levels(a2$model) == "model2"] <- "2"
levels(a2$model)[levels(a2$model) == "model2b"] <- "2 best"
levels(a2$model)[levels(a2$model) == "model3"] <- "3"
levels(a2$model)[levels(a2$model) == "model3b"] <- "3 best"
theme_update(strip.background = element_blank())
ggplot(a2) +
  geom_point(aes(Var1,  Var2, color = AVE_inner, alpha = AVE_inner)) +
  scale_color_viridis_c(begin = 0, end = 1) +
  labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
       size = "# genes") +
  guides(alpha = FALSE) +
  scale_alpha(range = c(0, 1)) +
  facet_wrap(~model)

a %>% 
  group_by(model) %>% 
  filter(AVE_inner == max(AVE_inner)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(Var1, Var2, size = AVE_inner, color = model, alpha = .1)) +
  labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
       size = "AVE_inner") +
  xlim(c(0, 1)) +
  ylim(c(0, 1))
a %>% 
  group_by(model) %>% 
  filter(AVE_inner == max(AVE_inner)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(Var1, Var2, size = AVE_inner, color = model, alpha = .1)) +
  labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
       size = "AVE_inner") +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  facet_wrap(~model)

model3b <- readRDS("model3_best.RDS")
a %>% 
  group_by(model) %>% 
  ungroup() %>% 
  filter(AVE_inner >= model3b$AVE$AVE_inner) %>% 
  ggplot() +
  geom_point(aes(Var1, Var2, size = AVE_inner, color = model)) +
  labs(x = "Shrinkage RNAseq", y = "Shrinkage 16S", 
       size = "AVE_inner") +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  facet_wrap(~model)
  

real <- dat %>% 
  filter(Var3 ==1)
r <- ggplot(real) +
  geom_point(aes(Var1, RNAseq, color = Var2, size = RNAseq)) +
  scale_color_viridis_c() +
  labs(x = "Shrinkage RNAseq", color = "Shrinkage 16S", title = "Effect of shrinkage parameter", 
       y = "# genes") +
  geom_vline(xintercept = 0.1)
s <- ggplot(real) +
  geom_point(aes(Var2, `16S`, color = Var1, size = RNAseq)) +
  scale_color_viridis_c() +
  labs(color = "Shrinkage RNAseq", x = "Shrinkage 16S", title = "Effect of shrinkage parameter", 
       y = "# OTUs")
shrin <- r + s




r2 <- ggplot(real) +
  geom_point(aes(RNAseq, AVE_inner,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_inner", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# genes", size = "Srhinkage 16S") +
  scale_size_continuous()
r3 <- ggplot(real) +
  geom_point(aes(`16S`, AVE_inner,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_inner", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# 16S", size = "Srhinkage 16S") +
  scale_size_continuous()
s2 <- ggplot(real) +
  geom_point(aes(RNAseq, AVE_outer,  color = Var1, size = Var2)) +
  scale_color_viridis_c() +
  labs(y = "AVE_outer", color = "Shrinkage RNAseq", title = "Effect of shrinkage parameter", 
       x = "# genes", size = "Srhinkage 16S") +
  scale_size_continuous()
AVEs <- r2 + s2

shrin/AVEs
