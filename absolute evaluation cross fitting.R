# Causal validation
# Absolute error
# TODO:
  # Absolute error can go negative. Take max{0, absolute error}? If we truncate the estimator at zero, will the CLT result still hold? Yes unless on the boundary?
  # nuisance.learner: method = "gradient boosting", "gradient boosting early stopping", "gradient boosting linear"
# For Weihan:
rm(list = ls())
source("~/Desktop/Research/Zijun/causal validation/helper.R")

n = 400  # number of observations
d = 5  # number of covariates

beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(1, d)
sigma = 1
n.fold = 5
method = "gradient boosting early stopping" # "linear", "gradient boosting early stopping"

# tau.hat function to evaluate
tau.hat.function = function(X){
  # set.seed(318)
  value = X %*% delta + delta0 # unbiased and noiseless
  value = X %*% (delta + rep(c(-0.5, 0.5), length = length(delta))) + delta0 # unbiased and noisy
  value = X %*% delta + delta0 + 2 # biased and noiseless
  value = X %*% (delta + rep(c(-0.5,0.5), length = length(delta))) + delta0 + 2 # biased and noisy
  return(value)
}


m = 400
record = list()
record$error.oracle = 
  record$error.semiEfficient = record$sd.semiEfficient =
  record$error.IF = record$sd.IF =
  record$error.baseline = rep(0, m)

set.seed(318)
for (i in 1 : m){
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, 0.5) # binary treatment
  mu0 = beta0 + X %*% beta 
  tau = X %*% delta + delta0
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
  # Estimation of nuisance functions
  # mu0.hat = mu0 + rnorm(n, 0, 5 / sqrt(n)) # unbiased is quite important
  # mu1.hat = mu1 + rnorm(n, 0, 5 / sqrt(n))
  # tau.hat = mu1.hat - mu0.hat
  # prop.hat = p
  folds = sample(rep(1 : n.fold, length.out = n))
  nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  record$error.oracle[i] = mean((tau.hat.function(X) - tau)^2)
  
  # Proposal based on the derived EIF
  IF.semiEfficient = tau.hat.function(X)^2 + second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, prop.hat = nuisance.hat$prop.hat, IF = T) - 2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, prop.hat = nuisance.hat$prop.hat, weight = tau.hat.function(X), IF = T)
  record$error.semiEfficient[i] = mean(IF.semiEfficient)
  record$sd.semiEfficient[i] = sd(IF.semiEfficient) / sqrt(n)
  
  # Method in Theorem 2 of "Validating Causal Inference Models via Influence Functions"
  B = 2 * W * (W - nuisance.hat$prop.hat) / nuisance.hat$prop.hat / (1 - nuisance.hat$prop.hat)
  A = W - nuisance.hat$prop.hat
  IF.IF = (1 - B) * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)^2 + B * Y * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X)) - A * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat -  tau.hat.function(X))^2 + tau.hat.function(X)^2
  record$error.IF[i] = mean(IF.IF)
  record$sd.IF[i] = sd(IF.IF) / sqrt(n)
  
  # Baseline method
  IF.baseline = (tau.hat.function(X) - (W * (Y - nuisance.hat$mu1.hat) / nuisance.hat$prop.hat + nuisance.hat$mu1.hat - (1 - W) * (Y - nuisance.hat$mu0.hat) / (1 - nuisance.hat$prop.hat) - nuisance.hat$mu0.hat))^2
  record$error.baseline[i] = mean(IF.baseline)
  
  if(i %*% 10 == 0){print(i)}
}

# results
data.frame(mean = unlist(lapply(record[-c(2,4)], mean)),
           sd = unlist(lapply(record[-c(2,4)], sd)),
           estimated.sd = c(NA, unlist(lapply(record[c(2,4)], mean)), NA))
