library("RGCCA")
library("ggplot2")
source("helper_functions.R")

otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)

A <- list(stools = otus_s, intestinal = otus_i)
C <- matrix(0, ncol = 2, nrow = 2, dimnames = list(names(A), names(A)))

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

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst, 
     # sgcca.factorial = sgcca.factorial)


samples <- data.frame(Stools = sgcca.horst$Y[[1]][, 1],
                 Intestinal = sgcca.horst$Y[[2]][, 1])

samples <- cbind(samples, meta)
subSamples <- samples[samples$Time == "T106", ]

ggplot(subSamples, aes(Stools, Intestinal)) +
  geom_text(aes(color =  Patient_ID, shape = Treatment, label = Sample_Code)) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("Samples integration", 
          subtitle = "Showing all samples after two years") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)")
  # coord_cartesian(ylim=c(-0.5, 0.5))
  # xlim(c(-0.25, 0.25)) + 
  # ylim(c(-0.5, 0.25))

variables <- data.frame(comp1 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 1]})),
                        comp2 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 2]})),
                        Origin = rep(names(A), sapply(A, ncol)))
variables$var <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
# Remove the variables that in both components are 0
keepComp1 <- abs(variables$comp1) > mean(abs(variables$comp1))
keepComp2 <- abs(variables$comp2) > mean(abs(variables$comp2))
subVariables <- variables[keepComp1 & keepComp2, ]

ggplot(subVariables, aes(comp1, comp2), color = Origin) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
  geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
  geom_text(aes(color = Origin, label = var)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  coord_cartesian(xlim=c(-0.25 , 0.25), ylim = c(-0.25, 0.25)) + 
  ggtitle("Variables important for the first two components", 
          subtitle = "Integrating stools and mucosa samples")

## Find the otus that are equivalent between datasets
library("dplyr")
comb <- expand.grid(rownames(otus_tax_i), rownames(otus_tax_s))
colnames(comb) <- c("intestinal", "stools")

eq <- apply(comb, 1, function(x){
  y <- x[2]
  x <- x[1]
  # If there is any NA then they are nor precise enough to say they are the same
  all(otus_tax_i[x, ] == otus_tax_s[y, ]) 
})

eqOTUS <- comb[eq & !is.na(eq), ]

subVariables2 <- subVariables
subVariables2$color <- "black"
library("plyr")
# mapvalues()

# ggplot(subVariables, aes(comp1, comp2), color = Origin) +
#   geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
#   geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
#   geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
#   geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
#   geom_text(aes(color = Origin, label = var)) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   coord_cartesian(xlim=c(-0.25 , 0.25), ylim = c(-0.25, 0.25)) + 
#   ggtitle("Variables important for the first two components", 
#           subtitle = "Integrating stools and mucosa samples")