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
alpha = 0.1 # 1 - alpha confidence level

# CI.lower and CI.upper is the lower and upper boundaries of the confidence interval of the comparison of LASSO V.S. Boosting. If CI.upper < 0, then LASSO is significantly better than Boosting.
selection.accuracy = function(CI.lower, CI.upper){
  selection = c("LASSO", "Equal", "Boosting")[1 + (CI.lower >= 0) + (CI.upper >= 0)]
  return(selection)
}

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
  plot.data = data.frame(selection.accuracy = c(unlist(lapply(record$record.absolute.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$plug.in[, 1] - x$plug.in[, 2] - qnorm(1 - alpha/2) * (x$sd.plug.in[, 1] + x$sd.plug.in[, 2]), CI.upper = x$plug.in[, 1] - x$plug.in[, 2] + qnorm(1 - alpha/2) * (x$sd.plug.in[, 1] + x$sd.plug.in[, 2])))})), 
                                   unlist(lapply(record$record.absolute.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$AVDS[, 1] - x$AVDS[, 2] - qnorm(1 - alpha/2) * (x$sd.AVDS[, 1] + x$sd.AVDS[, 2]), CI.upper = x$AVDS[, 1] - x$AVDS[, 2] + qnorm(1 - alpha/2) * (x$sd.AVDS[, 1] + x$sd.AVDS[, 2])))})),
                                   unlist(lapply(record$record.absolute.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$semiEfficient[, 1] - x$semiEfficient[, 2] - qnorm(1 - alpha/2) * (x$sd.semiEfficient[, 1] + x$sd.semiEfficient[, 2]), CI.upper = x$semiEfficient[, 1] - x$semiEfficient[, 2] + qnorm(1 - alpha/2) * (x$sd.semiEfficient[, 1] + x$sd.semiEfficient[, 2])))})),
                                   unlist(lapply(record$record.relative.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$semiEfficient - qnorm(1 - alpha/2) * x$sd.semiEfficient, CI.upper = x$semiEfficient + qnorm(1 - alpha/2) * x$sd.semiEfficient))})),
                                   unlist(lapply(record$record.relative.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$AVDS - qnorm(1 - alpha/2) * x$sd.AVDS, CI.upper = x$AVDS + qnorm(1 - alpha/2) * x$sd.AVDS))})),
                                   unlist(lapply(record$record.relative.total, function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$semiEfficient - qnorm(1 - alpha/2) * x$sd.semiEfficient, CI.upper = x$semiEfficient + qnorm(1 - alpha/2) * x$sd.semiEfficient))}))
                                   ))
  plot.data$Method = factor(rep(rep(c("plug in", "IF", "EIF"), times = rep(m, 3)), 2), levels = c("plug in", "IF", "EIF"))
  plot.data$Estimator = factor(rep(c("Absolute error", "Relative error"), times = rep(3 * m, 2)), levels = c("Absolute error", "Relative error"))

  line.width = 2; point.size = 3
  # LASSO
  # pdf(file = paste(plotDirectory, "/", plot.name, "_selection_accuracy", "_absolute_error", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[plot.data$Estimator == "Absolute error", ], aes(x = Method, y = selection.accuracy, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(title = "Absolute-error based",
         x = "Method", y = "Selection accuracy") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_x_discrete(labels = c("plug in", "IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  # dev.off()
  
  # Boosting
  # pdf(file = paste(plotDirectory, "/", plot.name, "_selection_accuracy", "_relative_error", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[(plot.data$Estimator == "Relative error") & (plot.data$Method != "plug in"), ], aes(x = Method, y = selection.accuracy, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(title = "Relative-error based",
         x = "Method", y = "Selection accuracy") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_x_discrete(labels = c("IF", "EIF/plug in")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  # dev.off()
}

