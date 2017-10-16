# Load libraries and files
library("RGCCA")
library("ggplot2")
source("helper_functions.R")


# Read files
otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)
tax_i <- read.csv(file = "intestinal_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)
tax_s <- read.csv(file = "stools_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)

# Prepare the metadata for the RGCCA package
 
# see https://stackoverflow.com/a/16200415/2886003
design <- model.matrix(~ . + 0, data=meta[, sapply(meta, is.factor)], 
             contrasts.arg = lapply(meta[, sapply(meta, is.factor)], contrasts, contrasts=FALSE))
stop("Continue here")
# Prepare input for the sgcca function
A <- list(stools = otus_s, intestinal = otus_i, metadata = meta)
# The design
C <- matrix(0, ncol = 2, nrow = 2, dimnames = list(names(A), names(A)))
C <- subSymm(C, "metadata", "stools", 1)
C <- subSymm(C, "metadata", "intestinal", 1)

# Shrinkage 
(shrinkage <- sapply(A, tau.estimate))

ncomp <- 2
ncomp <- rep(ncomp, length(A))

sgcca.centroid <-  sgcca(A, C, c1 = shrinkage,
                         ncomp = ncomp,
                         scheme = "centroid",
                         scale = TRUE,
                         verbose = FALSE)
names(sgcca.centroid$Y) <- names(A)
names(sgcca.centroid$a) <- names(A)
names(sgcca.centroid$astar) <- names(A)

sgcca.factorial <-  sgcca(A, C, c1 = shrinkage,
                          ncomp = ncomp,
                          scheme = "factorial",
                          scale = TRUE,
                          verbose = FALSE)
names(sgcca.factorial$Y) <- names(A)
names(sgcca.factorial$a) <- names(A)
names(sgcca.factorial$astar) <- names(A)

sgcca.horst <-  sgcca(A, C, c1 = shrinkage,
                      ncomp = ncomp,
                      scheme = "horst",
                      scale = TRUE,
                      verbose = FALSE)
names(sgcca.horst$Y) <- names(A)
names(sgcca.horst$a) <- names(A)
names(sgcca.horst$astar) <- names(A)

# Graphical explorations
theme_set(theme_bw())