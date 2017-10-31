# Helper functions and variables

library("ggplot2")
library("RGCCA")

#' Create a circle
#' @param center The position where the center of the circle is
#' @param diameter Arbitrary units of diameter of the circle
#' @param npoints Number of points of the circle (aka: definition of the circle)
#' @return A data.frame with the position of the points
circleFun <- function(center = c(-1,1), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circle <- circleFun(c(0, 0), 2, npoints = 100)

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
  colnames(M) <- rownames(M) = names(A)
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


#' Clean and prepare the data from IMNGS
#' 
#' Divides the taxonomy into a new matrix for each otu
#' 
#' @param taxonomy Last column of files a string ; separated with domain, 
#' phylum, vlass, order, family, genus and species.
#' @param otus The name of the rows
#' 
#' @return 
#' A matrix with the taxonomic information ready for the package phylo
taxonomy <- function(taxonomy, otus){
  taxonomy <- sapply(taxonomy, strsplit, split = ";")
  names(taxonomy) <- otus
  otus_tax <- t(sapply(taxonomy, '[', seq(max(sapply(taxonomy, length)))))
  colnames(otus_tax) <- c("Domain", "Phylum", "Class", "Order", 
                          "Family", "Genus", "Species")
  # Remove spaces
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "\\s", replacement = "")
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "[;:]", replacement = "")
  otus_tax <- apply(otus_tax, 1:2, sub, pattern = "^([a-z]__)", replacement = "")
  otus_tax[otus_tax == ""] <- NA # Remove empty cells
  otus_tax
}


# Check the taxonomy
# https://stackoverflow.com/q/7943695/2886003
#' Check if a vector is in the matrix
#' @param x The matrix to check if it is in the matrix
#' @param matrix The matrix where x is looked up.
#' 
#' @return A logical vector saying if the vector is in the matrix
fastercheck <- function(x, matrix){
  nc <- ncol(matrix)
  rec.check <- function(r, i, id){
    id[id] <- matrix[id, i] %in% r[i]
    if (i < nc & any(id)) 
      rec.check(r, i+1, id) 
    else 
      any(id)
  }
  apply(x, 1, rec.check, 1, rep(TRUE, nrow(matrix)))
}

tol21rainbow <-  c("#771155", "#AA4488", "#CC99BB", 
                   "#114477", "#4477AA", "#77AADD", 
                   "#117777", "#44AAAA", "#77CCCC", 
                   "#117744", "#44AA77", "#88CCAA", 
                   "#777711", "#AAAA44", "#DDDD77", 
                   "#774411", "#AA7744", "#DDAA77", 
                   "#771122", "#AA4455", "#DD7788")

colors <- c("#a692d2", "#6de14d", "#5a3bcb", "#b2e145", "#b844dd", "#50a93e",
            "#c640b6", "#6de697", "#472383", "#ddce3e", "#8766d6", "#8a9a37",
            "#517cd1", "#e24428", "#7dddcf", "#da425f", "#56ad7b", "#df4e98",
            "#c9e393", "#6a2d6d", "#de8d2f", "#3d3b76", "#d9ba73", "#d57ed2",
            "#3a682c", "#93295e", "#8f9c73", "#3c2032", "#d5d0bf", "#782e27",
            "#7ab4d5", "#ba5132", "#4d847e", "#de8f7d", "#2c3a29", "#d6abca",
            "#615325", "#ab6583", "#a47735", "#506280", "#89685f")

#' Calculates the angle between to slopes
#' 
#' @param x,y Slope of the lines
#' @note The default compares the first slope with the slope 1
#' @return The smaller angle between the slopes in degrees
angle <- function(x, y = 1){
  atan(abs((x-y)/(1+y*x)))*180/pi
}

#' Calculate the distance between a line and a point
#' 
#' The line is defined by points b and d.
#' @param p Point c(x, y)
#' @param b,d Points c(x, y) defining the line to calculate the distance with.
#' @return The units of distance between the point and the line
#' @note Change the d point to change the direction of the diagnoal
dist2d <- function(p, b = c(0, 0), d = c(1, 1)) {
  v1 <- b - d
  v2 <- p - b
  m <- cbind(v1, v2)
  abs(det(m))/sqrt(sum(v1*v1))
}

today <- format(Sys.time(), "%Y%m%d")

theme_set(theme_bw())
