# Helper functions and variables

library("ggplot2")
library("RGCCA")
library("STATegRa")
library("reshape2")

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


#' Compares the taxonomy of the otus
#' 
#' Given two taxonomy tables find which one is in which one
#' @param taxa_1, taxa_2 taxonomic table as described in taxonomy
#' @return A matrix of nrow(taxa_1) * nrow(taxa_2)
contingency_taxa <- function(taxa_1, taxa_2) {
  
  # Check if the input is factor or character
  txc1 <- apply(taxa_1, 2, is.character)
  txc2 <- apply(taxa_2, 2, is.character)
  txf1 <- apply(taxa_1, 2, is.factor)
  txf2 <- apply(taxa_2, 2, is.factor)
  
  if (!all(txc1 | txf1) | !all(txc2 | txf2)){
    stop("Taxa should be factors or characters (NAs allowed) ")
  }
  l <- sapply(rownames(taxa_2), function(x){
    nas_1 <- sum(is.na(taxa_1[x, ]))
    nas_2 <- sum(is.na(taxa_2[x, ]))
    (colSums(apply(taxa_1, 1,  `==`, taxa_2[x, ]), na.rm = TRUE) + 
        min(nas_2, nas_1))/7
  })
  l
}

#' Z-score to p-value
#' 
#' Calculates the p-value of a z-score
#' @param z the normalized value
#' @param one.sided Either NULL, - or +
#' @return The p-value 
convert.z.score <- function(z, one.sided = NULL) {
  # https://www.biostars.org/p/17227/#136907
  if (!is.null(one.sided) & !one.sided %in% c("+", "-")) {
    stop("one.sided should be NULL or + or -")
  }
  
  if (!is.numeric(z)) {
    stop("z should be numeric")
  }
  
  if(is.null(one.sided)) {
    pval <- pnorm(-abs(z))
    pval <- 2 * pval
  } else if(one.sided == "-") {
    pval <- pnorm(z)
  } else if (one.sided == "+") {
    pval <- pnorm(-z)
  }
  return(pval);
}

#' Calculates the z-score of two correlations
#' 
#' @param r1,r2 The correlation coefficients
#' @param n1,n2 The number of samples used to calculate the correlations
#' @return The z score of the comparison of the coefficients.
compare.correlations <- function(r1, r2, n1, n2) {
  r.norm <- function(r){log(abs((1+r)/(1-r)))/2}
  r1n <- r.norm(r1)
  r2n <- r.norm(r2)
  (r1n-r2n)/sqrt(1/(n1-3)+1/(n2-3))
}

#' Logical vectors of meta data
#' 
#' Given the names of the columns of the data calculates the logical vectors of 
#' each subset of the data
#' @param columns names of the columns to be used
#' @param data.frame with factors
#' @return a matrix with the logical values of each combination of the levels of 
#' the columns given for the data in each column. 
#' @note If some rows are all FALSE it means some values are NA.
allComb <- function(data, columns){
  if (is.null(dim(data))) {
    stop("data should be a data.frame or a matrix")
  }
  
  if (!is.character(columns)) {
    stop("columns should be a character")
  }
  
  if (!length(columns) >= 2) {
    stop("Several columns should be used")
  }
  
  keep <- columns %in% colnames(data)
  if (sum(keep) == 0) {
    stop("Names of columns not present on data")
  } else if (sum(keep) != length(keep)) {
    warning("Columns:", paste(columns[!keep]), "are not present on data")
  }
  
  data <- data[, columns[keep]]
  
  lvl <- sapply(data, levels)
  
  comb <- expand.grid(sapply(lvl, as.factor))
  comb2 <- apply(comb, 1, paste0, collapse = "_|_")
  out <- apply(comb, 1, function(x, data) {
    # Repeat the terms as much as the data
    combT <- sapply(x, function(y){rep(y, nrow(data))})
    # Compare the data with the levels
    check <- rowSums(data == combT)
    # Convert to logical by performing an &
    o <- check == 2
    o[is.na(o)] <- FALSE
    o
  }, data = data)
  colnames(out) <- comb2
  out
}
