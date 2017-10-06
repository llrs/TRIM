library("RGCCA")
library("ggplo2")
source("helper_functions.R")

otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")

A <- list(stools = otus_s, intestinal = otus_i)
C <- matrix(0, ncol = 2, nrow = 2, dimnames = list(names(A), names(A)))

subSymm <- function(m, x, y, val){
  m[x, y] <- val
  m[y, x] <- val
  m
}
C <- subSymm(C, "intestinal", "stools", 1)
# C <- subSymm(C, "metadata", "stools", 1)
# C <- subSymm(C, "metadata", "intestinal", 1)

(shrinkage <- sapply(A, tau.estimate))

ncomp <- c(2, 2)

sgcca.centroid <-  sgcca(A, C, c1 = shrinkage,
                     ncomp = ncomp,
                     scheme = "centroid",
                     scale = TRUE,
                     verbose = FALSE)

sgcca.factorial <-  sgcca(A, C, c1 = shrinkage,
                         ncomp = ncomp,
                         scheme = "factorial",
                         scale = TRUE,
                         verbose = FALSE)

sgcca.horst <-  sgcca(A, C, c1 = shrinkage,
                         ncomp = ncomp,
                         scheme = "horst",
                         scale = TRUE,
                         verbose = FALSE)

list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst, 
     sgcca.factorial = sgcca.factorial)