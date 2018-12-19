library("ggforce")
library("metagenomeSeq")
library("integration")
library("fgsea")
cd <- setwd("..")

# Load the helper file
today <- format(Sys.time(), "%Y%m%d")


intestinal <- "intestinal_16S"
rna <- "intestinal_RNAseq"

# Read the intestinal otus table
otus_table_i <- read.csv(
  file.path(intestinal, "OTUs-Table-new-biopsies.csv"),
  stringsAsFactors = FALSE, row.names = 1,
  check.names = FALSE
)
tax_i <- otus_table_i[, ncol(otus_table_i)]
otus_table_i <- otus_table_i[, -ncol(otus_table_i)]

# Extract the taxonomy and format it properly
otus_tax_i <- taxonomy(tax_i, rownames(otus_table_i))

# Load the input data
expr <- read.delim(file.path(rna, "taula_sencera2.tsv"), check.names = FALSE)

# Read the metadata for each type of sample
file_meta_i <- "intestinal_16S/db_biopsies_trim_seq16S_noBCN.txt"
meta_i <- read.delim(
  file_meta_i, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
file_meta_r <- file.path(rna, "metadata_25042018.csv")
meta_r <- read.delim(
  file_meta_r, check.names = FALSE,
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "")
)

setwd(cd)

# Correct the swapped samples and match metadata
expr <- norm_expr_colnames(expr)
# normalize names of samples
colnames(otus_table_i) <- gsub("[0-9]+\\.(.+)$", "\\1", colnames(otus_table_i))

# Normalize the RNA metadata
meta_r <- meta_r_norm(meta_r)

# Normalize the 16S intestinal metadata
meta_i <- meta_i_norm(meta_i)

# Check metadata with the names present in both datas
meta_r <- meta_r[meta_r$Seq_code_uDNA %in% colnames(otus_table_i) &
                   meta_r$`Sample Name_RNA` %in% colnames(expr), ]

# Subset the sequencing data
expr <- expr[, meta_r$`Sample Name_RNA`]
otus_table_i <- otus_table_i[, meta_r$Seq_code_uDNA]

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")
expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes=TRUE, log = TRUE)

# Filter expression
expr <- norm_RNAseq(expr_norm)

# Normalize OTUS
MR_i <- newMRexperiment(
  otus_table_i, 
  featureData = AnnotatedDataFrame(as.data.frame(otus_tax_i[rownames(otus_table_i), ]))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus_table_i <- MRcounts(MR_i, norm = TRUE, log = TRUE)

# Subset if all the rows are 0 and if sd is 0
otus_table_i <- otus_table_i[apply(otus_table_i, 1, sd) != 0, ]

# Create the different metadb matrices
Invariable <- cbind(mixOmics::unmap(meta_r$ID), 
                    Sex = mixOmics::unmap(meta_r$SEX)[, 1])
Location <- mixOmics::unmap(meta_r$Exact_location)
Location[is.na(Location)] <- 1
Time <- cbind(meta_r[, c("AGE_SAMPLE", "AgeDiag")],
              Transplant = mixOmics::unmap(meta_r$Transplant)[, 1],
              Surgery = mixOmics::unmap(meta_r$Surgery)[, 1])
Time$AgeDiag[is.na(Time$AgeDiag)] <- Time$AGE_SAMPLE[is.na(Time$AgeDiag)]

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus_table_i), 
          "Location" = Location, 
          "Time" = Time, 
          "Invariable" = Invariable)
A <- clean_unvariable(A)
# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "Time", 1)
C <- subSymm(C, "16S", "Location", 1)
C <- subSymm(C, "16S", "Invariable", 1)
C <- subSymm(C, "RNAseq", "Time", 1)
C <- subSymm(C, "RNAseq", "Location", 1)
C <- subSymm(C, "RNAseq", "Invariable", 1)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.249488046688595, 0, 1, 1, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage)
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

sgcca.centroid2 <- sgcca(
  A, D, c1 = shrinkage,
  ncomp = ncomp,
  scheme = "centroid",
  scale = TRUE,
  verbose = FALSE
)
sgcca.centroid2 <- improve.sgcca(sgcca.centroid2, names(A))
sgcca.centroid2$AVE

l2 <- lapply(A, function(x){
  cmdscale(dist(scale2(x)), k = 1)[, 1, drop = TRUE]
})

D <- C
for (i in seq_len(nrow(C))) {
  for (j in seq_len(i)) {
    if (i == j) {
      D[i, i] <- 0
    } else {
      
      D[i, j] <- cor(l2[[i]], l2[[j]])
      D[j, i] <- D[i, j]
    }
  }
}
D <- abs(D)
diag(D) <- 0
D <- subSymm(D, "Time", "Invariable", 1)
D <- subSymm(D, "Location", "Invariable", 1)
D <- subSymm(D, "Location", "Time", 0)
