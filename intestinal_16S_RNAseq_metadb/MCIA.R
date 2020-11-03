library("omicade4")
library("ggplot2")

# Load data ####

A <- readRDS("intestinal_16S_RNAseq_metadb/model3_TRIM.RDS")
meta_r <- readRDS( "intestinal_16S_RNAseq_metadb/meta.RDS")

l <- lapply(A[1:2], t)
# Check with mcia ####
out <- mcia(l)
# plot(out)
# plot(out$mcoa)
out$mcoa

# Extract results to compare with inteRmodels/RGCCA enhanced. ####
meta_r$Ileum <- ifelse(meta_r$Exact_location == "ILEUM", "ileum", "colon")
# From https://github.com/ComputationalSystemsBiology/momix-notebook/blob/master/scripts/runfactorization.R
data_plot <- cbind(out$mcoa$SynVar, meta_r) 
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = Ileum)) +
  theme_bw()
