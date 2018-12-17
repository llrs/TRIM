library("RGCCA")
library("gliomaData")

# Load data ####
data(ge_cgh_locIGR, package = "gliomaData")
A <- ge_cgh_locIGR$multiblocks
Loc <- factor(ge_cgh_locIGR$y)
levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
tau = c(1, 1, 0)

# rgcca algorithm using the dual formulation for X1 and X2 
# and the dual formulation for X3
A[[3]] = A[[3]][, -3]
C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)

# Compare features selected by different weights but same design ####

cmds <- lapply(A, function(x){cmdscale(dist(scale2(x)), k = 1)})

D2 <- D <- C
for (i in seq_len(nrow(C))) {
  Ai <- scale(A[[i]], scale = FALSE)
  for (j in seq_len(i)) {
    Bj <- scale(A[[j]], scale = FALSE)
    if (i == j) {
      D2[i, i] <- 0
    } else {
      D2[i, j] <- MatrixCorrelation::RV(Ai, Bj)
      D2[j, i] <- D2[i, j]
    }
  }
}

for (i in seq_len(nrow(C))) {
  for (j in seq_len(i)) {
    if (i == j) {
      D[i, i] <- 0
    } else {
      
      D[i, j] <- cor(cmds[[i]], cmds[[j]])
      D[j, i] <- D[i, j]
    }
  }
}
D <- abs(D)


result.sgcca = sgcca(A, D, c1 = c(.071,.2, 1), ncomp = c(2, 2, 2), 
                     scheme = "centroid", verbose = FALSE)
result.sgcca2 = sgcca(A, D2, c1 = c(.071,.2, 1), ncomp = c(2, 2, 2), 
                     scheme = "centroid", verbose = FALSE)
result.sgcca3 = sgcca(A, C, c1 = c(.071,.2, 1), ncomp = c(2, 2, 2), 
                     scheme = "centroid", verbose = FALSE)
result.sgcca4 = sgcca(A, D, c1 = c(.071,.2, 1), ncomp = c(2, 2, 2), 
                     scheme = "centroid", verbose = FALSE, init = "random")
comp.weight <- function(x, y){
  xdiff <- x != 0
  ydiff <- y != 0
  intersect(names(x[xdiff]), names(y[ydiff]))
}

l <- list(result.sgcca, result.sgcca2, result.sgcca3, result.sgcca4)
outm <- matrix(0, ncol = length(l), nrow = length(l))
for (i in seq_len(length(l))) {
  for (j in seq_len(i)) {
    if (i == j) {
      outm[i, i] <- 1
    } else {
      
      outm[i, j] <- length(comp.weight(l[[i]]$a[[1]][, 1], l[[j]]$a[[1]][, 1]))/145
      outm[j, i] <- outm[i, j]
    }
  }
}
outm

# Compare different weight but same design  ####

# Create the 1000 design matrix combinations
weight <- seq(from = 0, to = 1, by = .1)
C_list <- vector("list", length(weight))
names(C_list) <- as.character(weight)

for(i1 in weight){
  C_list[[as.character(i1)]] <- vector("list", length(weight))
  names(C_list[[as.character(i1)]]) <- as.character(weight)
  for (i2 in weight){
    C_list[[as.character(i1)]][[as.character(i2)]] <- vector("list", length(weight))
    names(C_list[[as.character(i1)]][[as.character(i2)]]) <- as.character(weight)
    for (i3 in weight) {
      C[2, 3] <- i1
      C[3, 2] <- i1
      C[1, 3] <- i2
      C[3, 1] <- i2
      C[1, 2] <- i3
      C[2, 1] <- i3
      C_list[[as.character(i1)]][[as.character(i2)]][[as.character(i3)]] <- C
    }
  }
}
C_list2 <- unlist(unlist(C_list, recursive = FALSE), recursive = FALSE)

# sgcca algorithm for each design
testing <- function(x) {
  result.sgcca <- sgcca(A, x, c1 = c(.071,.2, 1), ncomp = c(1, 1, 1), 
                       scheme = "centroid", verbose = FALSE)
  vs12 <- cor(result.sgcca$Y[[1]], result.sgcca$Y[[2]])
  vs13 <- cor(result.sgcca$Y[[1]], result.sgcca$Y[[3]])
  vs23 <- cor(result.sgcca$Y[[2]], result.sgcca$Y[[3]])
  cc <- x[1, 2]*abs(vs12)+x[1, 3]*abs(vs13) + x[2, 3]*abs(vs23)
  
  c(vs12 = vs12, vs13 = vs13, vs23 = vs23, 
    unlist(result.sgcca$AVE[c("AVE_inner", "AVE_outer")]), cc1 = cc)
}

out <- sapply(C_list2, testing)
out2 <- t(out)
var1 <- rep(weight, each = length(weight)^2, times = 1)
var2 <- rep(weight, each = length(weight), times = length(weight))
var3 <- rep(weight, each = 1, times = length(weight)^2)

# Design matrix:
# 0.0 var1 var2
# 0.0 0.0  var3
# 0.0 0.0  0.0
# load("../intestinal_16S_pathways_metadb/design_optimization.RData", verbose = TRUE)
def <- cbind.data.frame(out2, var1, var2, var3)
p <- ggplot(def)
p1 <- p + geom_point(aes(cc1, AVE_inner))
p2 <- p + geom_point(aes(cc1, AVE_outer))
p3 <- p + geom_count(aes(AVE_inner, AVE_outer))
q1 <- p + 
  geom_point(aes(var1, AVE_outer, col = "outer")) + 
  geom_smooth(aes(var1, AVE_outer, col = "outer")) +
  ylab("AVE")
q2 <- p + 
  geom_point(aes(var2, AVE_outer, col = "outer")) + 
  geom_smooth(aes(var2, AVE_outer, col = "outer")) + ylab("")

q3 <- p + 
  geom_point(aes(var3, AVE_outer)) + 
  geom_point(aes(var3, AVE_inner)) + ylab("") +
  geom_smooth(aes(var3, AVE_inner, col = "inner")) +
  geom_smooth(aes(var3, AVE_outer, col = "outer"))

r1 <- p +
  geom_smooth(aes(var1, AVE_inner, col = var2)) +
  geom_point(aes(var1, AVE_inner, col = var2))
r2 <- p +
  geom_smooth(aes(var1, AVE_inner, col = var3)) +
  geom_point(aes(var1, AVE_inner, col = var3))

r3 <- p +
  geom_smooth(aes(var2, AVE_inner, col = var1)) +
  geom_point(aes(var2, AVE_inner, col = var1))
r4 <- p +
  geom_smooth(aes(var2, AVE_inner, col = var3)) +
  geom_point(aes(var2, AVE_inner, col = var3))

r5 <- p +
  geom_smooth(aes(var3, AVE_inner, col = var2)) +
  geom_point(aes(var3, AVE_inner, col = var2))
r6 <- p +
  geom_smooth(aes(var3, AVE_inner, col = var1)) +
  geom_point(aes(var3, AVE_inner, col = var1))



library("patchwork")
(p1 + p2 + p3) / (q1 + q2 + q3 )

# Testing the tau effect on the AVE ####

min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
})
taus <- lapply(min_shrinkage, seq, to = 1, length.out = 10)
taus.combn <- expand.grid(taus)
out <- sapply(seq_len(nrow(taus.combn)), function(x){
  result.sgcca <- sgcca(A, C, c1 = taus.combn[x, ], ncomp = c(1, 1, 1), 
                        scheme = "centroid", verbose = FALSE)
  unlist(result.sgcca$AVE[c("AVE_inner", "AVE_outer")])
})
def <- cbind.data.frame(t(out), taus.combn)

# Add the default
result.sgcca = sgcca(A, C, c1 = c(.071,.2, 1), ncomp = c(1, 1, 1), 
                     scheme = "centroid", verbose = FALSE)
vec <- c(unlist(result.sgcca$AVE[c("AVE_inner", "AVE_outer")]), .071,.2, 1)
names(vec) <- colnames(def)
def <- rbind(def, vec)
p <- ggplot(def)
p1 <- p + geom_point(aes(GE, AVE_inner, col = "inner")) + ylab("AVE") +
  geom_smooth(aes(GE, AVE_inner, col = "inner")) +
  geom_smooth(aes(GE, AVE_outer, col = "outer")) +
  geom_point(aes(GE, AVE_outer, col = "outer")) + guides(col = FALSE)
p2 <- p + geom_point(aes(CGH, AVE_inner, col = "inner")) + ylab("") +
  geom_smooth(aes(CGH, AVE_inner, col = "inner")) +
  geom_smooth(aes(CGH, AVE_outer, col = "outer")) +
  geom_point(aes(CGH, AVE_outer, col = "outer")) + guides(col = FALSE)
p3 <- p + geom_point(aes(y, AVE_inner, col = "inner")) + ylab("") +
  geom_smooth(aes(y, AVE_inner, col = "inner")) +
  geom_smooth(aes(y, AVE_outer, col = "outer")) +
  geom_point(aes(y, AVE_outer, col = "outer")) + labs(col = 'AVE type') 

q <- ggplot(def) + geom_count(aes(AVE_outer, AVE_inner)) + 
  xlab("AVE outer") + ylab("AVE inner") + 
  scale_size_continuous(range = c(0.5, 2))

(p1 + p2 + p3) / (q)
