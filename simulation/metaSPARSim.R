library("metaSPARSim")

# Prepare data ####
A <- readRDS("intestinal_16S_RNAseq_metadb/model3_TRIM.RDS")
meta_r <- readRDS( "intestinal_16S_RNAseq_metadb/meta.RDS")
raw <- readRDS("intestinal_16S_RNAseq_integration/raw_otus.RDS")

# Estimate parameters
o <- estimate_parameter_from_data(raw_data = raw,
                                  norm_data = t(A[[2]]),
                                  conditions = split(seq_along(meta_r$ID), meta_r$ID))
# Error
names(o) <- names(split(seq_along(meta_r$ID), meta_r$ID))
sim_data <- metaSPARSim(o)

