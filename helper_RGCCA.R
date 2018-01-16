# Functions to use with RGCCA methods
library("ggplot2")

# I can't test if this is the McKeon homogeneity measure (no paper/reference)
# This shouldn't be calculated with a metablock
#' Calculates McKeon Homeogenity
#' 
#' Looks how much each blocks is related to the others blocks. This shouldn't 
#' be used with the metablock data (it is not a real block)
#' 
#' @param B The list of blocks
#' @param C The design matrix 
#' @return A matrix of the same dimensions as the list and design matrix with 
#' the "correlations" between blocks. 
#' @references There are no references other than the vignette of RGCCA
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
                 # It affects quite a lot between
                 # correlation c = 0 to covariation c = 1
                 c1 = rep(1, length(B)), 
                 ncomp = c(rep(2, (length(B)-1)),1),
                 # It doesn't seem to affect much (just the sign)
                 scheme = "horst", 
                 scale = TRUE,
                 verbose = FALSE)
  
  # Not sure of this simplification
  A <- sapply(seq_along(B), function(x){
    B[[x]][, unique(which(sgcca$a[[x]] != 0, arr.ind = TRUE)[, 1])]
  })
  
  J <- length(A)
  M <- matrix(0, J, J)
  rownames(M) <-  names(A)
  colnames(M) <- names(A)
  for (i in 1:J){
    for (j in seq_len(i)){
      M[i, j] <- rI(A[c(i, j)])
      M[j, i] <- M[i, j]
    }
  }
  dimnames(M) <- list(names(B), names(B))
  M
}

#' Subsitute in a symmetric matrix
#' 
#' @param m The symmetric matrix
#' @param x Row position
#' @param y Column position
#' @param val Value to insert in the given position
#' @return The symmetric matrix with the value inserted in the right positions
subSymm <- function(m, x, y, val){
  if (!isSymmetric(m)) {
    stop("m should be a symmetric matrix.")
  }
  m[x, y] <- val
  m[y, x] <- val
  m
}

#' Bootstrap sgcca
#' 
#' Performs the centroid bootstrap
#' 
#' @param A The list with the original data
#' @param C The symmetric matrix with the relationships between datsets.
#' @param shrinkage Shrinkage estimated (use the estimated for the original datastet)
#' @param nb_boot Number of bootstraps to perform
#' @return The outer weight of each variable of the input datasets. 
boot_sgcca <- function(A, C, shrinkage, nb_boot = 1000){
  STAB <- list()
  B <- A
  
  for (j in seq_along(A)) {
    STAB[[j]]<- matrix(NA, nb_boot, ncol(A[[j]]))
    colnames(STAB[[j]])<- colnames(A[[j]])
  }
  names(STAB) <- names(A)
  
  # Bootstrap the data
  for (i in 1:nb_boot){
    ind  <- sample(nrow(A[[1]]), replace = TRUE)
    Bscr <- lapply(A, `[`, i = ind)
    try( # Prevent the error from LAPACK subroutine
      {res <- sgcca(Bscr, C, c1 = shrinkage, 
                    ncomp = c(rep(1, length(A))),
                    scheme = "centroid", 
                    scale = TRUE)
      
      for (j in 1:J) {
        STAB[[j]][i, ] <- res$a[[j]]
      }}, silent = TRUE)
  }
  STAB
}

#' Evaluates the boostrapping of RGCCA
#' 
#' 
#' @param STAB List of weights of \code{rgcca} or \code{sgcca}
#' @return Lateral effect: Prints plots
boot_evaluate <- function(STAB){
  # Calculate how many are selected
  count <- lapply(STAB, function(x) {
    apply(x, 2, function(y){
      sum(y != 0, na.rm = TRUE)/(nb_boot - sum(is.na(STAB[[1]][, 1])))
    })
  })
  
  # Calculate the sign when selected
  sign <- lapply(STAB, function(x){colSums(sign(x), na.rm = TRUE)})
  
  # Calculate the mean and the standard error for each variable
  colMeAbs <- sapply(STAB, function(x){colMeans(abs(x), na.rm = TRUE)})
  seAbs <- sapply(STAB, function(x){
    apply(abs(x), 2, sd, na.rm = TRUE)/sqrt(nrow(x))
  })
  names(seAbs) <- names(STAB)
  names(colMeAbs) <- names(STAB)
  
  # Calculate the mean and the standard error for each variable
  colMe <- sapply(STAB, function(x){colMeans(x, na.rm = TRUE)})
  se <- sapply(STAB, function(x){
    apply(x, 2, sd, na.rm = TRUE)/sqrt(nrow(x))
  })
  names(se) <- names(STAB)
  names(colMe) <- names(STAB)
  
  # Merge the information in a table for each dataset
  var_info <- list(count, sign, colMeAbs, seAbs, colMe, se)
  consensus <- list()
  for (i in seq_along(STAB)){
    consensus[[i]] <- simplify2array(list("freq" = count[[i]], 
                                          "sign" = sign[[i]], 
                                          "colMeAbs" = colMeAbs[[i]], 
                                          "seAbs" = seAbs[[i]], 
                                          "colMe" = colMe[[i]], 
                                          "se" = se[[i]]))
    consensus[[i]] <- as.data.frame(consensus[[i]])
  }
  names(consensus) <- names(STAB)
  
  # Plot the summary of the bootstrapping
  for (i in seq_len(length(STAB))){
    p <- ggplot(consensus[[i]]) +
      geom_point(aes(sign, freq, col = colMeAbs, size = -log10(seAbs))) +
      ggtitle(paste("Selecting variable for", names(consensus)[i]))
    print(p)
    p <- ggplot(consensus[[i]]) +
      geom_point(aes(sign, freq, col = colMe, size = -log10(se))) +
      ggtitle(paste("Selecting variable for", names(consensus)[i]))
    print(p)
  }
}

