# plot for estimated absolute error of LASSO estimators with various hyper-parameters.
rm(list = ls())

path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")

library("ggplot2")
library("reshape2")

setting = c("similar HTE estimator")
plot.name = c("similar_HTE_estimator")

temp = readRDS(file.path(path, paste(setting, "absolute.rds", sep = " ")))
result = temp$result
lambda.seq = temp$lambda.seq

# Combine data into a data frame for ggplot
plot.data <- data.frame(
  lambda = lambda.seq,
  oracle.infData = result$oracle.infData,
  oracle = result$oracle,
  plug.in = result$plug.in,
  semiEfficient = result$semiEfficient,
  AVDS = result$AVDS,
  sd.semiEfficient = result$sd.semiEfficient,
  sd.AVDS = result$sd.AVDS,
  sd.plug.in = result$sd.plug.in
)

# Reshape into a long format for ggplot
plot.data.long <- melt(plot.data, id.vars = "lambda", 
                       measure.vars = c("oracle.infData", "oracle", "plug.in", "semiEfficient", "AVDS"),
                       variable.name = "Method", 
                       value.name = "Value")

# Plot
# pdf(file = paste(plotDirectory, "/", plot.name, "_absolute_error_curve", ".pdf", sep = ""), width = 4, height = 3.5)
g <- ggplot() +
  geom_line(data = subset(plot.data.long, Method == "oracle"), 
            aes(x = lambda, y = Value, color = "oracle"), 
            size = 1, linetype = "solid", show.legend = TRUE) +
  geom_line(data = subset(plot.data.long, Method == "semiEfficient"), 
            aes(x = lambda, y = Value, color = "EIF"), 
            size = 1, show.legend = FALSE) +
  geom_ribbon(data = subset(plot.data.long, Method == "semiEfficient"), 
              aes(x = lambda, ymin = Value - plot.data$sd.semiEfficient, 
                  ymax = Value + plot.data$sd.semiEfficient), 
              alpha = 0.2, color = NA,
              fill = "orange", show.legend = FALSE) +
  geom_line(data = subset(plot.data.long, Method == "AVDS"), 
            aes(x = lambda, y = Value, color = "IF"), 
            size = 1, show.legend = FALSE) +
  geom_ribbon(data = subset(plot.data.long, Method == "AVDS"), 
              aes(x = lambda, ymin = Value - plot.data$sd.AVDS, 
                  ymax = Value + plot.data$sd.AVDS), 
              alpha = 0.2, color = NA, 
              fill = "dark green", show.legend = FALSE) +
  geom_line(data = subset(plot.data.long, Method == "plug.in"), 
            aes(x = lambda, y = Value, color = "plug in"), 
            size = 1, show.legend = FALSE) +
  geom_ribbon(data = subset(plot.data.long, Method == "plug.in"), 
              aes(x = lambda, ymin = Value - plot.data$sd.plug.in, 
                  ymax = Value + plot.data$sd.plug.in), 
              alpha = 0.2, color = NA,
              fill = "dark blue", show.legend = FALSE) +
  labs(
    title = "", 
    x = expression(lambda), 
    y = "Estimated Absolute Error", 
    color = "Method"  # Title for the color legend
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("oracle" = "black", "plug in" = "dark blue", "EIF" = "orange", "IF" = "dark green")
  ) +
  theme(
    legend.position.inside = c(0.5, 5),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Print the plot
print(g)
# dev.off()



