path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")

library(ggplot2)

setting.total = "nuisance learner"  # "sample size", "nuisance learner" 
if(setting.total == "sample size"){
  setting.seq = paste("simulation", c(500, 1000, 1500, 2000, 2500, 3000))
}else if(setting.total == "nuisance learner"){
  setting.seq = paste("simulation", c("true", "linear", "ridge", "LASSO", "gradient boosting", "gradient boosting S learner"))  
}
record = readRDS(file.path(path, paste(setting.total, "rds", sep = ".")))
m = 100
alpha = 0.1 # 1 - alpha confidence level

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
  plot.data = data.frame(selection.accuracy = c(unlist(lapply(record$record.absolute.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$plug.in[, 1] - x$plug.in[, 2] - qnorm(1 - alpha/2) * (x$sd.plug.in[, 1] + x$sd.plug.in[, 2]), CI.upper = x$plug.in[, 1] - x$plug.in[, 2] + qnorm(1 - alpha/2) * (x$sd.plug.in[, 1] + x$sd.plug.in[, 2])))})), 
                                                unlist(lapply(record$record.absolute.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$AVDS[, 1] - x$AVDS[, 2] - qnorm(1 - alpha/2) * (x$sd.AVDS[, 1] + x$sd.AVDS[, 2]), CI.upper = x$AVDS[, 1] - x$AVDS[, 2] + qnorm(1 - alpha/2) * (x$sd.AVDS[, 1] + x$sd.AVDS[, 2])))})),
                                                unlist(lapply(record$record.absolute.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1, 1] - x$oracle[1, 2], CI.upper = x$oracle[1, 1] - x$oracle[1, 2]) == selection.accuracy(CI.lower = x$semiEfficient[, 1] - x$semiEfficient[, 2] - qnorm(1 - alpha/2) * (x$sd.semiEfficient[, 1] + x$sd.semiEfficient[, 2]), CI.upper = x$semiEfficient[, 1] - x$semiEfficient[, 2] + qnorm(1 - alpha/2) * (x$sd.semiEfficient[, 1] + x$sd.semiEfficient[, 2])))})),
                                                unlist(lapply(record$record.relative.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$semiEfficient - qnorm(1 - alpha/2) * x$sd.semiEfficient, CI.upper = x$semiEfficient + qnorm(1 - alpha/2) * x$sd.semiEfficient))})),
                                                unlist(lapply(record$record.relative.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$AVDS - qnorm(1 - alpha/2) * x$sd.AVDS, CI.upper = x$AVDS + qnorm(1 - alpha/2) * x$sd.AVDS))})),
                                                unlist(lapply(record$record.relative.total[[i]], function(x){mean(selection.accuracy(CI.lower = x$oracle[1], CI.upper = x$oracle[1]) == selection.accuracy(CI.lower = x$semiEfficient - qnorm(1 - alpha/2) * x$sd.semiEfficient, CI.upper = x$semiEfficient + qnorm(1 - alpha/2) * x$sd.semiEfficient))}))
  ))
  plot.data$Method = factor(rep(rep(c("plug in", "IF", "EIF"), times = rep(m, 3)), 2), levels = c("plug in", "IF", "EIF"))
  plot.data$Estimator = factor(rep(c("Absolute error", "Relative error"), times = rep(3 * m, 2)), levels = c("Absolute error", "Relative error"))
  
  line.width = 2; point.size = 3
  # LASSO
  pdf(file = paste(plotDirectory, "/", setting, " selection accuracy", " absolute error", ".pdf", sep = ""), width = 5, height = 3.5)
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
  dev.off()
  
  # Boosting
  pdf(file = paste(plotDirectory, "/", setting, " selection accuracy", " relative error", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[(plot.data$Estimator == "Relative error") & (plot.data$Method != "plug in"), ], aes(x = Method, y = selection.accuracy, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(title = "Relative-error based",
         x = "Method", y = "Selection accuracy") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_x_discrete(labels = c("IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
}
