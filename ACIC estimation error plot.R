# plot for sample size
path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")

library(ggplot2)

setting.seq = c("ACIC linear propensity linear HTE",
                "ACIC linear propensity nonlinear HTE",
                "ACIC nonlinear propensity linear HTE",
                "ACIC nonlinear propensity nonlinear HTE")

plot.name.seq = c("ACIC_linear_propensity_linear_HTE",
                "ACIC_linear_propensity_nonlinear_HTE",
                "ACIC_nonlinear_propensity_linear_HTE",
                "ACIC_nonlinear_propensity_nonlinear_HTE")

m = 100

custom_theme = theme_minimal(base_size = 14) + theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
  legend.direction = "horizontal",
  legend.position = "none" # "right",
  # panel.grid = element_blank()
)

# plot of error of estimated errors  
for(i in 1:length(setting.seq)){
    setting = setting.seq[i]
    plot.name = plot.name.seq[i]
    record = readRDS(file.path(path, paste(setting, "rds", sep = ".")))
    plot.data = data.frame(absolute.error = c(unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$plug.in[,1] - x$oracle[,1]))})), 
                                              unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$AVDS[,1] - x$oracle[,1]))})),  
                                              unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$semiEfficient[,1] - x$oracle[,1]))})),
                                              unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$plug.in[,2] - x$oracle[,2]))})), 
                                              unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$AVDS[,2] - x$oracle[,2]))})),                                                unlist(lapply(record$record.absolute.total, function(x){mean(abs(x$semiEfficient[,2] - x$oracle[,2]))})),
                                              unlist(lapply(record$record.relative.total, function(x){mean(abs(x$semiEfficient - x$oracle))})),
                                              unlist(lapply(record$record.relative.total, function(x){mean(abs(x$AVDS - x$oracle))})),
                                              unlist(lapply(record$record.relative.total, function(x){mean(abs(x$semiEfficient - x$oracle))}))))
    plot.data$Method = factor(rep(rep(c("plug in", "IF", "EIF"), times = rep(m, 3)), 3), levels = c("plug in", "IF", "EIF"))
    plot.data$Estimator = factor(rep(c("LASSO", "Boosting", "LASSO V.S. Boosting"), times = rep(3 * m, 3)), levels = c("LASSO", "Boosting", "LASSO V.S. Boosting"))
    plot.data$Method.Estimator = interaction(plot.data$Estimator, plot.data$Method)
    
  line.width = 2; point.size = 3
  # LASSO
  pdf(file = paste(plotDirectory, "/", plot.name, "_estimator_error", "_LASSO", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[plot.data$Estimator == "LASSO", ], aes(x = Method, y = absolute.error, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "LASSO",
         x = "Method", y = "Error of estimated validation error") +
    coord_cartesian(ylim = c(0, 150)) + 
    scale_x_discrete(labels = c("plug in", "IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
  
  # Boosting
  pdf(file = paste(plotDirectory, "/", plot.name, "_estimator_error", "_Boosting", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[plot.data$Estimator == "Boosting", ], aes(x = Method, y = absolute.error, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Boosting",
         x = "Method", y = "Error of estimated validation error") +
    coord_cartesian(ylim = c(0, 150)) + 
    scale_x_discrete(labels = c("plug in", "IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
  
  
  # LASSO V.S. Boosting
  pdf(file = paste(plotDirectory, "/", plot.name, "_estimator_error", "_LASSO_V.S._Boosting", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[(plot.data$Estimator == "LASSO V.S. Boosting") & (plot.data$Method != "plug in"), ], aes(x = Method, y = absolute.error, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "LASSO V.S. Boosting",
         x = "Method", y = "Error of estimated validation error") +
    coord_cartesian(ylim = c(0, 150)) + 
    scale_x_discrete(labels = c("IF", "EIF/plug in")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
}

