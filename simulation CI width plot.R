path = "~/Desktop/Research/Zijun/causal validation/Causal-validation/data"
plotDirectory = file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation")

library(ggplot2)

setting.total = "nuisance learner"  # "sample size", "nuisance learner" 
if(setting.total == "sample size"){
  setting.seq = paste("simulation", c(500, 1000, 1500, 2000, 2500, 3000))
  plot.name.seq = paste("simulation", c(500, 1000, 1500, 2000, 2500, 3000), sep = "_")
}else if(setting.total == "nuisance learner"){
  setting.seq = paste("simulation", c("true", "linear", "ridge", "LASSO", "gradient boosting", "gradient boosting S learner"))  
  plot.name.seq = paste("simulation", c("true", "linear", "ridge", "LASSO", "gradient_boosting", "gradient_boosting_S_learner"), sep = "_")  
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
  plot.name = plot.name.seq[i]
  plot.data = data.frame(width = c(unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.plug.in[, 1])})), 
                                   unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.AVDS[, 1])})),  
                                   unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.semiEfficient[, 1])})),  
                                   unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.plug.in[, 2])})), 
                                   unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.AVDS[, 2])})), 
                                   unlist(lapply(record$record.absolute.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.semiEfficient[, 2])})),  
                                   unlist(lapply(record$record.relative.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.semiEfficient)})),  
                                   unlist(lapply(record$record.relative.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.AVDS)})),  
                                   unlist(lapply(record$record.relative.total[[i]], function(x){mean(qnorm(1 - alpha / 2) * x$sd.semiEfficient)}))))
  plot.data$Method = factor(rep(rep(c("plug in", "IF", "EIF"), times = rep(m, 3)), 3), levels = c("plug in", "IF", "EIF"))
  plot.data$Estimator = factor(rep(c("LASSO", "Boosting", "LASSO V.S. Boosting"), times = rep(3 * m, 3)), levels = c("LASSO", "Boosting", "LASSO V.S. Boosting"))
  plot.data$Method.Estimator = interaction(plot.data$Estimator, plot.data$Method)
  
  line.width = 2; point.size = 3
  # LASSO
  pdf(file = paste(plotDirectory, "/", plot.name, "_CI_width", "_LASSO", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[plot.data$Estimator == "LASSO", ], aes(x = Method, y = width, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "HTE estimator LASSO",
         x = "Method", y = "Width") +
    coord_cartesian(ylim = c(0, 3.5)) + 
    scale_x_discrete(labels = c("plug in", "IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
  
  # Boosting
  pdf(file = paste(plotDirectory, "/", plot.name, "_CI_width", "_Boosting", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[plot.data$Estimator == "Boosting", ], aes(x = Method, y = width, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "HTE estimator Boosting",
         x = "Method", y = "Width") +
    coord_cartesian(ylim = c(0, 3.5)) + 
    scale_x_discrete(labels = c("plug in", "IF", "EIF")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
  
  
  # LASSO V.S. Boosting
  pdf(file = paste(plotDirectory, "/", plot.name, "_CI_width", "_LASSO_V.S._Boosting", ".pdf", sep = ""), width = 5, height = 3.5)
  g = ggplot(plot.data[(plot.data$Estimator == "LASSO V.S. Boosting") & (plot.data$Method != "plug in"), ], aes(x = Method, y = width, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "LASSO V.S. Boosting",
         x = "Method", y = "Width") +
    coord_cartesian(ylim = c(0, 3.5)) + 
    scale_x_discrete(labels = c("IF", "EIF/plug in")) +
    scale_fill_manual(values = c("plug in" = "navy",  "IF" = "dark green", "EIF" = "coral")) +  
    theme_minimal() +
    custom_theme
  print(g)
  dev.off()
}

# plot of choosing the better model with significance (significance level: 90%)

