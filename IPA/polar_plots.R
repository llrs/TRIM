# library(plotrix)
# library(gplots)
library(ggplot2)
library(forcats)

# setwd("IPA")
samples <- list.files(pattern = "*.csv")

fct_num <- function(x){
  as.numeric(levels(x))[x]
}

circular_plot <- function(x) {
  x <- x[x$Keep == "YES", ]
  x <- droplevels(x)
  x$`z-score` <- fct_num(x$`z-score`)
  x$color <- as.factor(ifelse(x$`z-score` > 0, "positive", "negative"))

  
  if (is(x$`-log(p-value)`, "factor")) {
    x$`-log(p-value)` <- fct_num(x$`-log(p-value)`)
  } else {
    x$`-log(p-value)` <- x$`-log(p-value)`
  }
  x <- x[order(abs(x$`-log(p-value)`), decreasing = TRUE), ]
  x$`Ingenuity Canonical Pathways` <- factor(x$`Ingenuity Canonical Pathways`,
                                             levels = x$`Ingenuity Canonical Pathways`)
  if (ncol(x) >= 9) {
    colnames(x) <- c("IPA", "log10pvalue", "ratio", "zscore", "molecules", 
                     "pvalue", "zs", "Keep", "color")
  } else {
    colnames(x) <- c("IPA", "log10pvalue", "ratio", "zscore", "molecules", 
                     "pvalue", "Keep", "color")
  }
  
  x$group <- "group_x"
  ggplot(x) + 
    geom_line(aes(IPA, log10pvalue, group = group), col = "grey") +
    geom_point(aes(IPA, log10pvalue, size = abs(zscore), col = color)) +
    coord_polar(theta = "x", direction = 1) +
    labs(x = "", y = "", color = "Sign", size = "z-score") +
    # ylim(c(0, 3)) +
    theme_minimal()  %+replace%
    theme(
      axis.text.y = element_text(angle = 0),
      axis.text = element_text(margin = margin(t = 20, b = 10), size = 7),
      axis.text.x = element_blank()
    ) +
    scale_alpha(range = c(0.5, 1))
  
  
}

sapply(samples, function(x) {
  name <- gsub("canonical_(.+)\\.csv", "\\1", x)
  y <- read.csv(x, check.names = FALSE)
  ggsave(filename = paste0(name, "_white.png"), circular_plot(y))
})


