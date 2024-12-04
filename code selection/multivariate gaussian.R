# Load necessary library
# type I error
library(MASS)

# Set the seed for reproducibility
set.seed(123)

d = 40
p.value = rep(0, d)
# Specify the mean vector for the multivariate normal
mu <- rep(0, d); mu[1] = -0.1 # Four-dimensional mean

# Specify the covariance matrix for the correlated entries
Sigma = matrix(0.8, nrow = d, ncol = d)
diag(Sigma) = 1

# Generate multivariate Gaussian data
n <- 100  # Number of samples

m = 100
record = list(); record$p.value = record$min.gap = rep(0, m)
for(i in 1:m){
  data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  # Convert to data frame for better viewing
  data_df <- as.data.frame(data)
  names(data_df) <- paste("X", seq(1, d), sep = "")
  
  # Display the first few rows of the data
  test.stats = apply(data_df, 2, mean)
  winner = which.min(test.stats) # the smaller, the better
  
  # p values
  # naive
  var.hat = apply(data_df - data_df[, rep(winner, d)], 2, sd) / sqrt(n)
  p.value = pnorm((test.stats[winner] - test.stats[-winner]) / sqrt(var.hat[-winner]), 0, 1)
  
  record$p.value[i] = max(p.value)
  # condition on nuisance  
  record$min.gap[i] = max((test.stats[winner] - test.stats[-winner]) / sqrt(var.hat[-winner]))
  
}

par(mfrow = c(2, 2))
# quite large, but <= 0.5
hist(record$p.value, xlim = c(0, 1))
hist(record$min.gap)


