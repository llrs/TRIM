
corr <- read.csv("20181119_genus_correlation_all_model2_best.csv")


# Find if it has too many good correlations
# Meaning that it has an outlier
library("dplyr")
corr2 <- corr %>% 
  select(Microorganism, Gene, Correlation, pvalue) %>% 
  group_by(Microorganism) %>% 
  count() %>% 
  filter(n > 60) %>% 
  ungroup()

corr2 <- corr %>% 
  filter(Microorganism %in% corr2$Microorganism) %>% 
  droplevels()
dfs <- split(corr2, corr2$Microorganism)
library("ggplot2")
theme_set(theme_bw())
corr2 %>% 
  ggplot() +
  geom_violin(aes(Microorganism, abs(Correlation)))
out <- sapply(dfs, function(x) {
  cors <- abs(x$Correlation)
  sum(cors > 0.4)/length(cors)
})
