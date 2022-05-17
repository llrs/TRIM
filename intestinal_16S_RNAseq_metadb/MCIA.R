library("omicade4")
library("ggplot2")

# Load data ####

A <- readRDS("model3_TRIM.RDS")
meta_r <- readRDS( "meta.RDS")

l <- lapply(A[1:2], t)
# Check with mcia ####
out <- mcia(l)
plot(out)
# plot(out$mcoa)
# out$mcoa

# Extract results to compare with inteRmodels/RGCCA enhanced. ####
meta_r$Ileum <- ifelse(meta_r$Exact_location == "ILEUM", "ileum", "colon")
# From https://github.com/ComputationalSystemsBiology/momix-notebook/blob/master/scripts/runfactorization.R
data_plot <- cbind(out$mcoa$SynVar, meta_r) 
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = Ileum, shape = IBD), size = 5) +
  theme_bw() +
  labs(title = "MCIA on HSCT")
ggsave("Figures/MCIA_SynVar.png")
ggsave(filename = "~/Documents/projects/thesis/images/hsct-mcia.png", width = 170,
       units = "mm", dpi = 300, bg = "white")
ggplot(data_plot) +
  geom_point(aes(SynVar1, SynVar2, col = IBD, shape = IBD)) +
  theme_bw()

library("pROC")
roc <- multiclass.roc(data_plot$SynVar1[!is.na(meta_r$Ileum)], response = meta_r$Ileum[!is.na(meta_r$Ileum)], levels = unique(meta_r$Ileum[!is.na(meta_r$Ileum)]))
auc(roc)

# So it seems that the synthetic scores is the one that has the common space
# But all the variables are all important for all the dimensions, making it hard to understant the meaning. 
vars <- data.frame(out$mcoa$Tco) %>% 
  mutate(variable = rownames(.),
         type = ifelse(startsWith(variable, "ENSG"), "gene", "micro"))
ggplot(vars) +
  geom_density(aes(SV1, col = type)) +
  theme_minimal() +
  labs(title = "Weights of variables for SynVar",
       subtitle = "All variables have a weight different from 0")
ggsave("Figures/MCIA_Tco.png")
ggplot(vars) +
  geom_density(aes(SV2, col = type))

# How do we compare RGCCA and MCIA? 
# We can't compare the genes and microorganisms selected because it selects everything.
# We could place a threshold but where would we set it? 
# What would be inner AVE:
C <- matrix(rep(1, 4), ncol = 2) # What is the model ?
sum(C*out$mcoa$cov2/2)/(sum(C)/2) # Seems more realistic but I don't understand what it is.
# Or 
RGCCA:::ave_inner(C, out$mcoa$SynVar)

