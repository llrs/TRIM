# Helper functions

library("ggplot2")
library("RGCCA")

# Circle
circleFun <- function(center = c(-1,1),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circle <- circleFun(c(0, 0), 2, npoints = 100)

# I can't test if this is the McKeon homogeneity measure
# This shouldn't be calculated with a metablock
McKeonHomeogenity <- function(B, C) {
  if (!all(length(B) == ncol(C) & ncol(C) == nrow(C))) {
    stop("Number of blocks and design matrix are not coherent")
  }
  
  # helper function provided in the vignette
  rI <- function(A){
    J <- length(A)
    res <- rgcca(A, scale = TRUE, verbose = FALSE)
    Y <- Reduce("cbind", res$Y)
    rI <- 1/(J-1)*(cov2(rowSums(Y))/sum(apply(Y, 2, cov2))-1)
    rI
  }
  
  sgcca <- sgcca(B, C,
                 c1 = rep(1, length(B)), # It affects quite a lot from correlation c == 0 to covariation c == 1
                 ncomp = c(rep(2, (length(B)-1)),1),
                 scheme = "horst", # It doesn't seem to affect much (just the sign)
                 scale = TRUE,
                 verbose = FALSE)
  
  # Not sure of this simplification
  A = sapply(seq_len(length(A)), function(x){
    B[[x]][, unique(which(sgcca$a[[x]]!=0, arr.ind = TRUE)[, 1])]
  })
  
  J = length(A)
  M = matrix(0, J, J)
  colnames(M) = rownames(M) = names(A)
  for (i in 1:J){
    for (j in seq_len(i)){
      M[i, j] = rI(A[c(i, j)])
      M[j, i] = M[i, j]
    }
  }
  dimnames(M) <- list(names(B), names(B))
  M
}
