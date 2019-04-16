library("ggplot2")
library("gganimate")

theme_update(strip.background = element_blank())

out <- readRDS("models_summary.RDS")
out$Interaction <- as.numeric(grepl("i$", out$Model))
out <- out[out$Interaction != 1 & out$Component == "comp1", ]
basic <- ggplot(out, aes(GE, M, col = IBD, group = 'Sample Name_RNA')) +
  geom_point()
basic + geom_line(colour = 'grey')
anim <- basic +
  transition_states(Model, wrap = FALSE) +
  # view_follow(fixed_x = c(-0.35, 0.2), fixed_y = c(-0.6, 1)) +
  ggtitle('Model {closest_state}')
animate(anim)
basic2 <- ggplot(out, aes(GE, M, shape = IBD, group = 'Sample Name_RNA', color = SESCD_local)) +
  geom_point()
anim2 <- basic2 +
  transition_states(Model, wrap = FALSE) +
  # view_follow(fixed_x = c(-0.35, 0.2), fixed_y = c(-0.6, 1)) +
  ggtitle('Model {closest_state}')
animate(anim2)
