library("RGCCA")
library("ggplot2")
source("helper_functions.R")

otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)

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
theme_set(theme_bw())

(shrinkage <- sapply(A, tau.estimate))

ncomp <- c(2, 2)

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

list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst, 
     sgcca.factorial = sgcca.factorial)


df <- data.frame(Stools = sgcca.horst$Y[[1]][, 1],
                 Intestinal = sgcca.horst$Y[[2]][, 1])

df <- cbind(df, meta)
subDf <- df[df$Time == "T26", ]

ggplot(subDf, aes(Stools, Intestinal)) +
  geom_text(aes(color =  Patient_ID, shape = Treatment, label = Sample_Code)) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
  # coord_cartesian(ylim=c(-0.5, 0.5))
  # xlim(c(-0.25, 0.25)) + 
  # ylim(c(-0.5, 0.25))



variables <- data.frame(comp1 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 1]})),
                        comp2 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 2]})),
                        orig = rep(names(A), sapply(A, ncol)))
# Remove the variables that in both components are 0
keep <-  apply(variables[, c("comp1", "comp2")], 1, function(x){all(x != 0)})
subVariables <- variables[keep, ]

ggplot(subVariables, aes(comp1, comp2), color = orig) +
  # geom_path(aes(x, y), data = circle) + 
  geom_text(aes(color = orig, label = rownames(subVariables)))
  # coord_cartesian(xlim=c(-0.3, 0.3), ylim = c(-0.3, 0.3))
