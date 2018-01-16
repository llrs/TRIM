cd <- setwd("..")

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
setwd(cd)

# Prepare the metadata for the RGCCA package

meta2 <- meta
keepCol <- sapply(meta, is.factor)
postTreatment <- c("Birth_date", "Sample_Code", "Patient_ID", "HSCT_responder")
postTreatment["Transplant"] <- TRUE
keepCol[postTreatment] <- FALSE
for (col in names(keepCol)[keepCol]){
  levels(meta2[, col]) <- seq_along(levels(meta2[, col]))
}

# see https://stackoverflow.com/a/16200415/2886003
# nas <- getOption("na.action")
# options(na.action = "na.pass")
# design <- model.matrix(~ . + 0, data = meta[, keepCol, drop = FALSE], 
#                        contrasts.arg = lapply(meta[, keepCol, drop = FALSE], 
#                                               contrasts, contrasts = FALSE))
# options(na.action = nas)

# For those without information we say they are not in any group
# design[is.na(design)] <- 0 
# attributes(design) <- attributes(design)[-4]

# Set metadb with a sigle variable with several options
metadb <- meta2[, keepCol]
metadb <- apply(metadb, 1:2, as.numeric)
metadb[is.na(metadb)] <- 0
# Prepare input for the sgcca function
A <- list(stools = otus_s, intestinal = otus_i, metadata = metadb)
# The design
C <- matrix(0, ncol = length(A), nrow = length(A), 
            dimnames = list(names(A), names(A)))
C <- subSymm(C, "metadata", "stools", 1)
C <- subSymm(C, "metadata", "intestinal", 1)


# Shrinkage 
(shrinkage <- sapply(A, tau.estimate))
(min_shrinkage <- sapply(A, function(x){1/sqrt(ncol(x))}))
# Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# We want to keep the covariance of the metadata, hence forcing the 1:
shrinkage[length(shrinkage)] <- 1

ncomp <- 2
ncomp <- rep(ncomp, length(A))

sgcca.centroid <-  sgcca(A, C, c1 = c(1, 1,1),
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


McKeonHomeogenity(A, C)

samples <- data.frame(Stools = sgcca.centroid$Y[[1]][, 1],
                      Intestinal = sgcca.centroid$Y[[2]][, 1])

names(colors) <- unique(meta$ID)
samples <- cbind(samples, meta)
samples$Patient_ID <- as.factor(samples$Patient_ID)

pdf(paste0("Figures/", today, "_plots.pdf"))

# Labels of the samples
label <- strsplit(as.character(samples$Sample_Code), split = "_")
labels <- sapply(label, function(x){
  if (length(x) == 5){
    paste(x[5], sep = "_")
    # x[5]
  }
  else if (length(x) != 5) {
    paste(x[4], sep = "_")
    # x[4]
  }
})
samples <- cbind(samples, labels)
samples$Time <- factor(samples$Time, levels(samples$Time)[c(1, 5, 6, 3, 4, 2)])

for (p in seq_along(levels(samples$Time))){
  a <- ggplot(samples, aes(Stools, Intestinal)) +
    geom_text(aes(color =  Patient_ID, 
                  label = paste(ID, labels, sep = "_"))) + 
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) + 
    xlab("Stools (component 1)") +
    ylab("Mucosa (component 1)") +
    guides(col = guide_legend(title="Patient ID")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) + 
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}
ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color =  ID, 
                label = paste(Time, labels, sep = "_"))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title="Patient ID")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = colors)

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color = HSCT_responder , 
                label = paste(ID, labels, sep = "_"))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title = "Responders")) + 
  theme(plot.title = element_text(hjust = 0.5))


ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color = Endoscopic_Activity, 
                label = paste(ID, labels, sep = "_"))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title = "Endoscopic Activity")) + 
  theme(plot.title = element_text(hjust = 0.5))

ggplot(samples, aes(Stools, Intestinal)) +
  geom_text(aes(color = Time, 
                label = paste(ID, labels, sep = "_"))) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") + 
  xlab("Stools (component 1)") +
  ylab("Mucosa (component 1)") +
  guides(col = guide_legend(title = "Time")) + 
  theme(plot.title = element_text(hjust = 0.5))

variables <- data.frame(comp1 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 1]})),
                        comp2 = unlist(sapply(sgcca.centroid$a, 
                                              function(x){x[, 2]})),
                        Origin = rep(names(A), sapply(A, ncol)))
variables$var <- gsub("^.*\\.(.*)$", "\\1", rownames(variables))

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


# To calculate the conficence interval on selecting the variable
# this interval should reduce as we fit a better model/relationship
nb_boot <- 1000 # number of bootstrap samples
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

# Calculate how many are selected
count <- lapply(STAB, function(x) {
  apply(x, 2, function(y){
    sum(y != 0)/nb_boot
  })
})

# Calculate the sign when selected
sign <- lapply(STAB, function(x){colSums(sign(x))})

# Calculate the mean and the standard error for each variable
colMeAbs <- sapply(STAB, function(x){colMeans(abs(x))})
seAbs <- sapply(STAB, function(x){
  apply(abs(x), 2, sd)/sqrt(nrow(x))
})
names(seAbs) <- names(STAB)
names(colMeAbs) <- names(STAB)

# Calculate the mean and the standard error for each variable
colMe <- sapply(STAB, function(x){colMeans(x)})
se <- sapply(STAB, function(x){
  apply(x, 2, sd)/sqrt(nrow(x))
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
for (i in seq_len(2)){
  p <- ggplot(consensus[[i]]) +
    geom_point(aes(sign, freq, col = colMeAbs, size = -log10(seAbs))) +
    ggtitle(paste("Selecting variable for", names(consensus)[i])) +
    xlab("Sign") + 
    ylab("Frequency")
  print(p)
  p <- ggplot(consensus[[i]]) +
    geom_point(aes(sign, freq, col = colMe, size = -log10(se))) +
    ggtitle(paste("Selecting variable for", names(consensus)[i])) + 
    xlab("Sign") + 
    ylab("Frequency")
  print(p)
}

# Select the most important variables of the bootstrap
selectedVar <- sapply(consensus[1:2], selectVar)

# Find the organisms most important and shared between stools and intestinal
tax_s_s <- unique(tax_s[selectedVar[["stools"]], ])
tax_i_s <- unique(tax_i[selectedVar[["intestinal"]], ])

# Find the OTUs names
s_in_i <- fastercheck(tax_s_s, tax_i_s)
com <- tax_i_s[s_in_i, ]
i <- tax_i_s[!s_in_i,]
s <- tax_s_s[!fastercheck(tax_i_s, tax_s_s), ]

dev.off()

# Write the output files
write.csv(com, file = "important_common_microrg.csv", 
          row.names = FALSE, na = "")
write.csv(i, file = "important_intestinal_microrg.csv", 
          row.names = FALSE, na = "")
write.csv(s, file = "important_stools_microrg.csv", 
          row.names = FALSE, na = "")