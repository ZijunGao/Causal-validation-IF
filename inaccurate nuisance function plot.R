# plot for sample size
path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")


setting = c("inaccurate nuisance function estimator absolute", "inaccurate nuisance function estimator relative")
record.absolute = readRDS(file.path(path, paste(setting[1], "rds", sep = ".")))
record.relative = readRDS(file.path(path, paste(setting[2], "rds", sep = ".")))

plot.name = c("inaccurate_nuisance_function_estimator_absolute", "inaccurate_nuisance_function_estimator_relative")

library(ggplot2) 

index = 1
custom_theme = theme_minimal(base_size = 14) + theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
  legend.position = "none",
  panel.grid = element_blank()
)

pdf(file = paste(plotDirectory, "/", plot.name[1], "_LASSO", ".pdf", sep = ""), width = 3.8, height = 3.5)
plot.data.LASSO = data.frame("absolute_error" = c(record.absolute$plug.in[index, 1], record.absolute$AVDS[index, 1], record.absolute$semiEfficient[index, 1]), 
                             "absolute_error_sd" = c(record.absolute$sd.plug.in[index, 1], record.absolute$sd.AVDS[index, 1], record.absolute$sd.semiEfficient[index, 1]))
plot.data.LASSO$method = factor(c("plug in", "IF", "EIF"), levels = c("plug in", "IF", "EIF"))
line.width = 2; point.size = 3
# Plot for tau.hat = LASSO
ggplot(plot.data.LASSO, aes(x = method, y = absolute_error, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = absolute_error - qnorm(0.95) * absolute_error_sd, ymax = absolute_error + qnorm(0.95) * absolute_error_sd), width = 0.2) +
  geom_hline(yintercept = record.absolute$oracle.infData[index, 1], linetype = "dashed", color = "dark red") +
  labs(title = "LASSO HTE estimator", 
       x = "Method", 
       y = "Estimated absolute error") +
  scale_color_manual(values = c("navy", "dark green", "coral"), guide = "none") +
  ylim(-12, 12) + 
  custom_theme
dev.off()


pdf(file = paste(plotDirectory, "/", plot.name[1], "_xgboost", ".pdf", sep = ""), width = 3.8, height = 3.5)
plot.data.xgboost = data.frame("absolute_error" = c(record.absolute$plug.in[index, 2], record.absolute$AVDS[index, 2], record.absolute$semiEfficient[index, 2]), 
                             "absolute_error_sd" = c(record.absolute$sd.plug.in[index, 2], record.absolute$sd.AVDS[index, 2], record.absolute$sd.semiEfficient[index, 2]))
plot.data.xgboost$method = factor(c("plug in", "IF", "EIF"), levels = c("plug in", "IF", "EIF"))
line.width = 2; point.size = 3
# Plot for tau.hat = LASSO
ggplot(plot.data.xgboost, aes(x = method, y = absolute_error, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = absolute_error - qnorm(0.95) * absolute_error_sd, ymax = absolute_error + qnorm(0.95) * absolute_error_sd), width = 0.2) +
  geom_hline(yintercept = record.absolute$oracle.infData[index, 2], linetype = "dashed", color = "dark red") +
  labs(title = "Boosting HTE estimator", 
       x = "Method", 
       y = "Estimated absolute error") +
  scale_color_manual(values = c("navy", "dark green", "coral")) +
  ylim(-12, 12) + 
  custom_theme
dev.off()


pdf(file = paste(plotDirectory, "/", plot.name[2], ".pdf", sep = ""), width = 3.8, height = 3.5)

plot.data.relative = data.frame("relative_error" = c(-record.relative$AVDS[index], -record.relative$semiEfficient[index]), 
                               "relative_error_sd" = c(record.relative$sd.AVDS[index], record.relative$sd.semiEfficient[index]))
plot.data.relative$method = factor(c("IF", "EIF/plug in"), levels = c("IF", "EIF/plug in"))
line.width = 2; point.size = 3
# Plot for tau.hat = LASSO
ggplot(plot.data.relative, aes(x = method, y = relative_error, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = relative_error - qnorm(0.95) * relative_error_sd, ymax = relative_error + qnorm(0.95) * relative_error_sd), width = 0.2) +
  geom_hline(yintercept = -record.relative$oracle.infData[index], linetype = "dashed", color = "dark red") +
  labs(title = "LASSO V.S. Boosting", 
       x = "Method", 
       y = "Estimated relative error") +
  scale_color_manual(values = c("dark green", "coral")) + 
  ylim(-12, 12) + 
  custom_theme
dev.off()
