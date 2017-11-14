# Load the helper file
source("helper_functions.R")

# Read files
otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)
tax_i <- read.csv(file = "intestinal_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)
tax_s <- read.csv(file = "stools_16S/taxonomy.csv", 
                  row.names = 1, stringsAsFactors = FALSE)

# Remove outlier See PCA
keep <- !grepl("28_T52_T_DM_CH", meta$Sample_Code) 
meta <- meta[keep, ]
otus_s <- otus_s[keep, ]
otus_i <- otus_i[keep, ]

# Subset the ileum
keepIleum <- meta$CD_Aftected_area == "ILEUM"
meta <- meta[keepIleum, ]
otus_s <- otus_s[keepIleum, ]
otus_i <- otus_i[keepIleum, ]

# Remove empty otus
removeOtus_i <- colSums(otus_i) != 0
removeOtus_s <- colSums(otus_s) != 0

otus_s <- otus_s[, removeOtus_s]
otus_i <- otus_i[, removeOtus_i]

tax_s <- tax_s[removeOtus_s, ]
tax_i <- tax_i[removeOtus_i, ]

##### RGCCA #####
A <- list(stools = otus_s, intestinal = otus_i)
C <- matrix(0, ncol = length(A), nrow = length(A), 
            dimnames = list(names(A), names(A)))

# Move to its own folder
setwd("ileum_integration")

C <- subSymm(C, "intestinal", "stools", 1)

(shrinkage <- sapply(A, tau.estimate))
(min_shrinkage <- sapply(A, function(x){1/sqrt(ncol(x))}))
# Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

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


samples <- data.frame(Stools = sgcca.centroid$Y[[1]][, 1],
                      Intestinal = sgcca.centroid$Y[[2]][, 1])
lmt <- lm(Intestinal ~ Stools, data = samples)

if (lmt$coefficients[2] > 0) {
  d <- c(1, 1)
} else if (lmt$coefficients[2] < 0) {
  d <-  c(-1, -1)
}

dist <- apply(samples, 1, dist2d, d = d)

names(colors) <- unique(meta$Patient_ID)
samples <- cbind(samples, meta, "dist" = dist)
samples$Patient_ID <- as.factor(samples$Patient_ID)

pdf(paste0("Figures/", today, "_RGCCA_plots.pdf"))

# Plot if the coherence between samples has a specific pattern
ggplot(samples) + 
  geom_point(aes(Patient_ID, log10(dist), col = Involved_Healthy)) + 
  facet_grid(~ Time)

# Labels of the samples
label <- strsplit(as.character(samples$Sample_Code), split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[1], x[2], x[5], sep = "_")
    # x[5]
  }
  else if (length(x) != 5) {
    paste(x[1], x[4], sep = "_")
    # x[4]
  }
})

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color =  Patient_ID, label = labels)) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle(paste0("Samples by time")) + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors) + 
  geom_abline(intercept = 0, slope = d[1], linetype = 2) + 
  facet_grid(~Time, scales = "free")

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color =  Patient_ID, label = labels)) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle(paste0("Samples by time")) + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors) + 
  geom_abline(intercept = 0, slope = d[1], linetype = 2) + 
  facet_grid(~Time)

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color =  Patient_ID, label = labels)) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors) + 
  geom_abline(intercept = 0, slope = d[1], linetype = 2)

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
comb <- expand.grid(rownames(tax_i), rownames(tax_s), stringsAsFactors = FALSE)
colnames(comb) <- c("intestinal", "stools")

eq <- apply(comb, 1, function(z){
  y <- z[2]
  x <- z[1]
  # If there is any NA then they are nor precise enough to say they are the same
  # OTU
  all(tax_i[x, ] == tax_s[y, ]) 
})

# Create a matrix with TRUE or FALSE when two otus are equivalent
dim(eq) <- c(nrow(tax_i), nrow(tax_s))
eq[is.na(eq)] <- FALSE

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time", "Active_area"))
subCors <- sapply(as.data.frame(time_area), function(x){
  if (sum(x) == 0) {
    return()
  }
  i <- otus_i[x, ]
  s <- otus_s[x, ]
  subCor <- cor(i, s, use = "pairwise.complete.obs", 
                method = "spearman")
  y <- subCor[eq]
  y[!is.na(y)]
}, simplify = FALSE)

subCors <- subCors[colSums(time_area) >= 4]
corEqOtus <- melt(subCors)
name <- strsplit(corEqOtus$L1, "_|_")
corEqOtus$Time <- sapply(name, function(x){x[1]})
corEqOtus$`Active area` <- sapply(name, function(x){x[3]})
# Set the factors of Time in the order we want
corEqOtus$Time <- factor(corEqOtus$Time, 
                         levels = c("T0", "TM36", "TM48", "T26", "T52", "T106"))

ggplot(corEqOtus) +
  geom_boxplot(aes(`Active area`, value)) + 
  facet_grid(~Time) +
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools")


# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x){x[, 1]})
# Select the most important variables
selectedVar <- sapply(comp1, function(x){
  q <- quantile(x, na.rm = TRUE)
  if (q["25%"] < q["75%"] ) {
    names(x)[(x < q["25%"] | x > q["75%"]) &  !is.na(x)]
  } else if (q["25%"] > q["75%"] ) {
    names(x)[(x > q["25%"] | x < q["75%"]) &  !is.na(x)]
  } else if (q["25%"] == q["75%"] ) {
    names(x)[x != 0]
  }
})

# Find the organisms most important and shared between stools and intestinal
tax_s_s <- tax_s[selectedVar[["stools"]], ]
tax_i_s <- tax_i[selectedVar[["intestinal"]], ]
s_in_i <- fastercheck(unique(tax_s_s), unique(tax_i_s))
com <- unique(tax_i_s)[s_in_i, ]
i <- unique(tax_i_s[!s_in_i,])
s <- unique(tax_s_s)[!fastercheck(unique(tax_i_s), unique(tax_s_s)), ]
write.csv(com, file = "important_common_microrg.csv", 
          row.names = FALSE, na = "")
write.csv(i, file = "important_intestinal_microrg.csv", 
          row.names = FALSE, na = "")
write.csv(s, file = "important_stools_microrg.csv", 
          row.names = FALSE, na = "")

comp1 <- sapply(comp1, '[', seq(max(sapply(comp1, length))))
rownames(comp1) <- seq_len(nrow(comp1))

# Plot the densities of the loadings
comp1 <- melt(comp1)[2:3]
colnames(comp1) <- c("Origin", "Loadings")
ggplot(comp1) +
  geom_density(aes(x = Loadings, 
                   y = ..scaled..,
                   fill = Origin), alpha = 0.5) +
  ggtitle("Importance of the otus of each data set") +
  ylab("Scaled density") +
  xlab("OTUs weight") +
  facet_grid(~Origin) + 
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x){x[, 2]})
comp2 <- sapply(comp2, '[', seq(max(sapply(comp2, length))))
rownames(comp2) <- seq_len(nrow(comp2))

comp2 <- melt(comp2)[2:3]
colnames(comp2) <- c("Origin", "Loadings")
ggplot(comp2) +
  geom_density(aes(x = Loadings, y = ..scaled.., fill = Origin), alpha = 0.5) +
  ggtitle("Importance of each block variable", 
          subtitle = "Second component") +
  ylab("Scaled density") +
  xlab("OTUs weight") +
  facet_grid(~Origin) + 
  guides(fill = FALSE) 

# To calculate the conficence interval on selecting the variable
# this interval should reduce as we fit a better model/relationship
nb_boot <- max(ncol(otus_i), ncol(otus_s)) # number of bootstrap samples
J <- length(A)
STAB <- list()
B <- lapply(A, cbind)

for (j in 1:J) {
  STAB[[j]]<- matrix(0, nb_boot, NCOL(A[[j]]))
  colnames(STAB[[j]])<- colnames(B[[j]])
}
names(STAB) <- names(B)

# Bootstrap the data
for (i in 1:nb_boot){
  ind  <- sample(NROW(B[[1]]), replace = TRUE)
  Bscr <- lapply(B, function(x) x[ind, ])
  res <- sgcca(Bscr, C, c1 = shrinkage, 
               ncomp = c(rep(1, length(B))),
               scheme = "centroid", 
               scale = TRUE)
  
  for (j in 1:J) {
    STAB[[j]][i, ] <- res$a[[j]]
  }
}

# Calculate the mean and the standard error for each variable
colMe <- sapply(STAB, colMeans)
se <- sapply(STAB, function(x){
  apply(x, 2, sd)/sqrt(nrow(x))
})
names(se) <- names(STAB)
names(colMe) <- names(STAB)
# Select the block we want to plot the variables for
for (i in seq_along(se)) {
  a <- cbind("SE" = se[[i]], "mean" = colMe[[i]])
  a <- as.data.frame(a)
  a <- a[order(a$mean, a$SE, decreasing = c(TRUE, FALSE)), ]
  
  p <- ggplot(a) + 
    geom_pointrange(aes(x = 1:nrow(a), y = mean, 
                        ymin = mean - SE, ymax = mean + SE)) + 
    ggtitle(names(se)[i]) + 
    xlab("Features") +
    ylab("Weight")
  print(p)
}

dev.off()

##### STATegRa #####

keepT0 <- meta$Time == "T106"
keep <- keepT0

# Create ExpressionSet objects
expr <- as.matrix(t(otus_i))[, keep]
colnames(expr) <- rownames(meta[keep, ])

eS_i <- createOmicsExpressionSet(expr, pData = meta[keep, ])
# eS_i <- ExpressionSet(assayData = as.matrix(otus_i))
expr <- as.matrix(t(otus_s))[, keep]
colnames(expr) <- rownames(meta[keep, ])

eS_s <- createOmicsExpressionSet(expr, pData = meta[keep, ])
# eS_s <- ExpressionSet(assayData = as.matrix(otus_s))

pdf(paste0("Figures/", today, "ILEUM_STATegRa_plots.pdf"))

# Selecting components
cc <- selectCommonComps(t(otus_i[keep, ]), t(otus_s[keep, ]), Rmax = 3)
PCA.selection(t(otus_i[keep, ]), fac.sel = "single%", varthreshold = 0.03)$numComps
PCA.selection(t(otus_s[keep, ]), fac.sel = "single%", varthreshold = 0.03)$numComps
(ms <- modelSelection(list(eS_i, eS_s), Rmax = 7, fac.sel = "single%",
                      varthreshold = 0.03))

plot(cc$pssq)
plot(cc$pratios)
# Omics Integration
discoRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
                              Names = c("Intestinal", "Stools"),
                              method = "DISCOSCA",
                              Rcommon = ms$common,
                              Rspecific = ms$dist,
                              center = TRUE, scale = TRUE)
plotVAF(discoRes)

jiveRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
                             Names = c("Intestinal", "Stools"),
                             method = "JIVE",
                             Rcommon = ms$common,
                             Rspecific = ms$dist,
                             center=TRUE, scale=TRUE)
o2plsRes <- omicsCompAnalysis(list("Intestinal" = eS_i, "Stools" = eS_s), 
                              Names = c("Intestinal", "Stools"),
                              method="O2PLS", 
                              Rcommon = ms$common,
                              Rspecific = ms$dist,
                              center=TRUE, scale=TRUE, weight=TRUE)
dev.off() 
save.image(file = paste0(today, "_ileum_integration.RData"))
