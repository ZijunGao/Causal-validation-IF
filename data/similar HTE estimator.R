# plot for sample size
path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")


setting = c("similar HTE estimator absolute", "similar HTE estimator relative")
record.absolute = readRDS(file.path(path, paste(setting[1], "rds", sep = ".")))
record.relative = readRDS(file.path(path, paste(setting[2], "rds", sep = ".")))



library(ggplot2)

custom_theme = theme_minimal(base_size = 14) + theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
  legend.direction = "vertical",
  legend.position = "right",
  panel.grid = element_blank()
)

index = 26 # 26
lambda.index = 5; lambda.index.reference = 8
pdf(file = paste(plotDirectory, "/similar HTE estimator", " LASSO", ".pdf", sep = ""), width = 5, height = 3.5)
plot.data = data.frame(
  Method = factor(rep(c("plug in", "IF", "EIF"), each = 3), levels = c("plug in", "IF", "EIF")),
  Estimator= factor(rep(c("HTE estimator 1", "HTE estimator 2", "HTE estimator 1 V.S. HTE estimator 2"), times = 3), 
                    levels = c("HTE estimator 1", "HTE estimator 2", "HTE estimator 1 V.S. HTE estimator 2")),
  error = c(record.absolute$plug.in[index, lambda.index],
            record.absolute$plug.in[index, lambda.index.reference],
            record.relative$semiEfficient[index, lambda.index],
            record.absolute$AVDS[index, lambda.index],
            record.absolute$AVDS[index, lambda.index.reference],
            record.relative$AVDS[index, lambda.index],
            record.absolute$semiEfficient[index, lambda.index],
            record.absolute$semiEfficient[index, lambda.index.reference],
            record.relative$semiEfficient[index, lambda.index]), 
  sd = c(record.absolute$sd.plug.in[index, lambda.index],
         record.absolute$sd.plug.in[index, lambda.index.reference],
         record.relative$sd.semiEfficient[index, lambda.index],
         record.absolute$sd.AVDS[index, lambda.index],
         record.absolute$sd.AVDS[index, lambda.index.reference],
         record.relative$sd.AVDS[index, lambda.index],
         record.absolute$sd.semiEfficient[index, lambda.index],
         record.absolute$sd.semiEfficient[index, lambda.index.reference],
         record.relative$sd.semiEfficient[index, lambda.index])
)
line.width = 2; point.size = 3

# Plot for tau.hat = LASSO

# Plot the data using ggplot2
ggplot(plot.data[plot.data$Method != "plug in",], aes(x = Method, y = pmax(0, error), color = Method, shape = Estimator)) +  
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = pmax(0, error - qnorm(0.9) * sd), ymax = pmax(0, error + qnorm(0.9) * sd)), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = record.relative$oracle.infData[index, lambda.index], linetype = "dashed", color = "dark red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "", y = "Error evaluation", x = "Method") +
  scale_shape_manual(values = c(16, 17, 18), 
                     labels = c(latex2exp::TeX("$\\hat{\\tau}_1$"), latex2exp::TeX("$\\hat{\\tau}_1$"), latex2exp::TeX("$\\hat{\\tau}_1$ V.S. $\\hat{\\tau}_2$"))) +
  scale_color_manual(values = c("dark green", "coral")) +
  # ylim(-2.5, 5.1) + 
  custom_theme
# dev.off()

