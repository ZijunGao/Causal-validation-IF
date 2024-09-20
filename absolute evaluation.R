# causal validation
# absolute evaluation
# TODO:
source("~/Desktop/Research/Zijun/causal validation/helper.R")

n = 100  # number of observations
d = 5  # number of covariates

beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(1, d)
sigma = 1

# predicted function
tau.hat.function = function(X){
  value = X %*% delta + delta0 + rnorm(dim(X)[1], 0, 1 / sqrt(n)) # unbiased and noiseless
  value = X %*% delta + delta0 + rnorm(dim(X)[1], 0, 2) # unbiased and noisy
  value = X %*% delta + delta0 + 2 + rnorm(dim(X)[1], 0, 1 / sqrt(n)) # biased and noiseless
  value = X %*% delta + delta0 + 2 + rnorm(dim(X)[1], 0, 2) # biased and noisy
  return(value)
}

m = 100
record = list()
record$error.oracle = 
  record$error.semiEfficient = 
  record$error.baseline = rep(0, m)
set.seed(318)
for (i in 1 : m){
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, 0.5)  # binary treatment
  mu0 = beta0 + X %*% beta 
  tau = X %*% delta + delta0
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
  # nuisance estimation for validation
  mu0.hat = mu0 + rnorm(n, 0, 5 / sqrt(n)) # unbiased is quite important
  mu1.hat = mu1 + rnorm(n, 0, 5 / sqrt(n))
  tau.hat = mu1.hat - mu0.hat
  prop.hat = p
  
  record$error.oracle[i] = mean((tau.hat.function(X) - tau)^2)
  record$error.semiEfficient[i] = mean((tau.hat.function(X))^2) + second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = mu1.hat, mu0.hat = mu0.hat, prop.hat = prop.hat) - 2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = mu1, mu0.hat = mu0, prop.hat = prop.hat, weight = tau.hat.function(X))
  record$error.baseline[i] = mean((tau.hat.function(X) - (W * (Y - mu1.hat) / prop.hat + mu1.hat - (1 - W) * (Y - mu0.hat) / (1 - prop.hat) - mu0.hat))^2)
}

# results
data.frame(mean = unlist(lapply(record, mean)),
           sd = unlist(lapply(record, sd)))

