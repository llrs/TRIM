cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")
library("integration")

# Read files
otus_i <- read.csv(file = "intestinal_16S/otus_coherent.csv")
otus_s <- read.csv(file = "stools_16S/otus_coherent.csv")
meta <- read.csv(file = "meta_coherent.csv", row.names = 1)
tax_i <- read.csv(
  file = "intestinal_16S/taxonomy.csv",
  row.names = 1, stringsAsFactors = FALSE
)
tax_s <- read.csv(
  file = "stools_16S/taxonomy.csv",
  row.names = 1, stringsAsFactors = FALSE
)
eqOTUS <- read.csv("equivalent_otus.csv", stringsAsFactors = FALSE)
setwd(cd)
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
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)

# Move to its own folder
C <- subSymm(C, "intestinal", "stools", 1)

(shrinkage <- sapply(A, tau.estimate))
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

sgcca.centroid <- sgcca(
  A, C, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid <- improve.sgcca(sgcca.centroid, names(A))

samples <- data.frame(
  Stools = sgcca.centroid$Y[[1]][, 1],
  Intestinal = sgcca.centroid$Y[[2]][, 1]
)
lmt <- lm(Intestinal ~ Stools, data = samples)

if (lmt$coefficients[2] > 0) {
  d <- c(1, 1)
} else if (lmt$coefficients[2] < 0) {
  d <- c(-1, -1)
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
labels <- sapply(label, function(x) {
  if (length(x) == 5) {
    paste(x[2], x[5], sep = "_")
    # x[5]
  }
  else if (length(x) != 5) {
    paste(x[4], sep = "_")
    # x[4]
  }
})

samples$Time <- factor(samples$Time, levels(samples$Time)[c(1, 5, 6, 3, 4, 2)])
for (p in seq_along(levels(samples$Time))) {
  a <- ggplot(samples, aes(Stools, Intestinal)) +
    geom_text(aes(color = Patient_ID, label = labels)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) +
    xlab("Stools (component 1)") +
    ylab("Mucosa (component 1)") +
    guides(col = guide_legend(title = "Patient ID")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    geom_abline(intercept = 0, slope = d[1], linetype = 2) +
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

for (p in seq_along(levels(samples$Time))) {
  a <- ggplot(samples, aes(Stools, Intestinal)) +
    geom_text(aes(color = Patient_ID, label = labels)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) +
    xlab("Stools (component 1)") +
    ylab("Mucosa (component 1)") +
    guides(col = guide_legend(title = "Patient ID")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    geom_abline(intercept = 0, slope = d[1], linetype = 2) +
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color = Patient_ID, label = labels)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title = "Patient ID")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors) +
  geom_abline(intercept = 0, slope = d[1], linetype = 2)

variables <- data.frame(
  comp1 = unlist(sapply(
    sgcca.centroid$a,
    function(x) {
      x[, 1]
    }
  )),
  comp2 = unlist(sapply(
    sgcca.centroid$a,
    function(x) {
      x[, 2]
    }
  )),
  Origin = rep(names(A), sapply(A, ncol))
)
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
  coord_cartesian(xlim = c(-0.25, 0.25), ylim = c(-0.25, 0.25)) +
  ggtitle(
    "Variables important for the first two components",
    subtitle = "Integrating stools and mucosa samples"
  )


# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x) {
  x[, 1]
})
cvariables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {
  x[, 2]
})
variables_weight(comp2)

# To calculate the conficence interval on selecting the variable
# this interval should reduce as we fit a better model/relationship
# Bootstrap of sgcca
boot <- boot_sgcca(A, C, shrinkage, 1000)

saveRDS(boot, file = "bootstrap.RDS")

# Evaluate the boostrap effect and plot
boot_evaluate(boot$STAB)

## Use the previously calculated table of equivalent OTUs
u_i <- unique(eqOTUS$intestinal)
u_s <- unique(eqOTUS$stools)

### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time"))
subCors <- sapply(as.data.frame(time_area), function(x) {
  cor(
    otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs",
    method = "spearman"
  )
}, simplify = FALSE)

# Tranform the data for the plot
cors <- sapply(subCors, function(x) {
  x <- x[upper.tri(x)]
  x[!is.na(x)]
})

cors <- cors[colSums(time_area) >= 4]
corEqOtus <- melt(cors)

# Set the factors of Time in the order we want
ggplot(corEqOtus) +
  geom_boxplot(aes(L1, value)) +
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools") +
  xlab("CD afected area")


### Compare the equivalent otus in different settings
time_area <- allComb(meta, c("Time", "Active_area"))
subCors <- sapply(as.data.frame(time_area), function(x) {
  cor(
    otus_i[x, u_i], otus_s[x, u_s], use = "pairwise.complete.obs",
    method = "spearman"
  )
}, simplify = FALSE)


# Tranform the data for the plot
cors <- sapply(subCors, function(x) {
  x <- x[upper.tri(x)]
  x[!is.na(x)]
})

cors <- cors[colSums(time_area) >= 4]
corEqOtus <- melt(cors)

name <- strsplit(corEqOtus$L1, "_|_")
corEqOtus$Time <- sapply(name, function(x) {
  x[1]
})
corEqOtus$Site <- sapply(name, function(x) {
  x[3]
})
# Set the factors of Time in the order we want
corEqOtus$Time <- factor(
  corEqOtus$Time,
  levels = c("T0", "TM36", "TM48", "T26", "T52", "T106")
)

ggplot(corEqOtus) +
  geom_boxplot(aes(Site, value)) +
  facet_grid(~Time) +
  ylab("Correlation with stools") +
  ggtitle("Comparison with stools")


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
(ms <- modelSelection(
  list(eS_i, eS_s), Rmax = 7, fac.sel = "single%",
  varthreshold = 0.03
))

plot(cc$pssq)
plot(cc$pratios)
# Omics Integration
discoRes <- omicsCompAnalysis(
  list("Intestinal" = eS_i, "Stools" = eS_s),
  Names = c("Intestinal", "Stools"),
  method = "DISCOSCA",
  Rcommon = ms$common,
  Rspecific = ms$dist,
  center = TRUE, scale = TRUE
)
plotVAF(discoRes)

jiveRes <- omicsCompAnalysis(
  list("Intestinal" = eS_i, "Stools" = eS_s),
  Names = c("Intestinal", "Stools"),
  method = "JIVE",
  Rcommon = ms$common,
  Rspecific = ms$dist,
  center = TRUE, scale = TRUE
)
o2plsRes <- omicsCompAnalysis(
  list("Intestinal" = eS_i, "Stools" = eS_s),
  Names = c("Intestinal", "Stools"),
  method = "O2PLS",
  Rcommon = ms$common,
  Rspecific = ms$dist,
  center = TRUE, scale = TRUE, weight = TRUE
)
dev.off()
