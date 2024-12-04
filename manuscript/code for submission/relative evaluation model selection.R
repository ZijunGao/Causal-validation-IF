# Causal Validation: Absolute Error for Model Selection
# Varying setup:
# p
# mu0
# tau
# Observations:
# green: sensitive to delta0
# TODO:
  # varying p
  # xgboost
  # relative performance: correlated CI
# Load required library
rm(list = ls())
library(glmnet)
source("helper.R")


# Apply the sequence to the same dataset and report errors
n = 400; # n = 400
n.train = 400
n.test = 10000; d = 5; n.fold = 2 # 5; 2
beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(0, d); delta[1:2] = rep(2, 2)
sigma = 2 # 1
p = 0.5 # 0.5; 0.2

lambda.seq = 2^(-seq(-1, 6) / 2) # simple: 2^(-seq(-1, 4) / 1); complex: 2^(-seq(-5, 20) / 5)
method = "gradient boosting S learner" # "linear"; "gradient boosting"; "gradient boosting early stopping"; "gradient boosting S learner" 
m = 100
record = list()
record$oracle.infData = record$oracle =  record$semiEfficient = record$sd.semiEfficient = record$AVDS = record$sd.AVDS = matrix(0 , nrow = m, ncol = length(lambda.seq) - 1)

setting = "similar HTE estimator" # "homophily", "similar HTE estimator"
if(setting == "homophily"){
  method = "linear"
  # use the quadratic HTE: tau = X %*% delta + delta0 + X[, 1]^2 
}
if(setting == "similar HTE estimator"){
  method = "linear"
  # use the quadratic HTE: tau = X %*% delta + delta0 + X[, 1]^2 
}

set.seed(318)
# Test data
X.test = matrix(rnorm(n.test * d), nrow = n.test)  # covariates
W.test = rbinom(n.test, 1, p) # binary treatment
mu0.test = beta0 + X.test %*% beta # + X.test[, 1]^2
tau.test = X.test %*% delta + delta0 + X.test[, 1]^2 # + X.test[, 1]^2 
mu1.test = mu0.test + tau.test
Y0.test = mu0.test + rnorm(n.test, 0, sigma)  # outcome with treatment effect
Y1.test = Y0.test + tau.test
Y.test = Y1.test * W.test + Y0.test * (1 - W.test)

for(i in 1 : m){
  # Generate  data
  # training
  X.train = matrix(rnorm(n.train * d), nrow = n.train)  # covariates
  W.train = rbinom(n.train, 1, p) # binary treatment
  mu0.train = beta0 + X.train %*% beta # + X.train[, 1]^2
  tau.train = X.train %*% delta + delta0 + X.train[, 1]^2 # + X.train[, 1]^2 
  mu1.train = mu0.train + tau.train
  Y0.train = mu0.train + rnorm(n.train, 0, sigma)  # outcome with treatment effect
  Y1.train = Y0.train + tau.train
  Y.train = Y1.train * W.train + Y0.train * (1 - W.train)
  
  # validation
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, p) # binary treatment
  mu0 = beta0 + X %*% beta # + X[, 1]^2
  tau = X %*% delta + delta0 + X[, 1]^2 # + X[, 1]^2 
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
  # Fit LASSO with the current lambda
  lasso.fit = glmnet(x = cbind(X.train, W.train, W.train * X.train), Y.train, alpha = 1, lambda = lambda.seq)
  
  # Define tau.hat.function based on the LASSO fit
  tau.hat.functions = function(X, j) {
    data.cnt = cbind(X, 0, 0 * X)
    data.trt = cbind(X, 1, 1 * X)
    tau.hat = predict(lasso.fit, newx = data.trt)[, j] - predict(lasso.fit, newx = data.cnt)[, j]
    return(tau.hat)
  }
  
  # MSE
  folds = sample(rep(1 : n.fold, length.out = n))
  nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  
  # Oracle
  record$oracle.infData[i, ] = sapply(seq(1, length(lambda.seq) - 1), function(x){mean((tau.hat.functions(X.test, x) - tau.test)^2)}) - mean((tau.hat.functions(X.test, length(lambda.seq)) - tau.test)^2)

  record$oracle[i, ] = sapply(seq(1, length(lambda.seq) - 1), function(x){mean((tau.hat.functions(X, x) - tau)^2)}) - mean((tau.hat.functions(X, length(lambda.seq)) - tau)^2)
  
  # Semi-parametrically efficient
  temp.semiEfficient = sapply(seq(1, length(lambda.seq) - 1), function(x){efficient.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, length(lambda.seq))})})
  record$semiEfficient[i,] = unlist(temp.semiEfficient[1,])
  record$sd.semiEfficient[i,] = unlist(temp.semiEfficient[2, ])
  
  # AVDS
  temp.AVDS = sapply(seq(1, length(lambda.seq) - 1), function(x){AVDS.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, length(lambda.seq))})})
  record$AVDS[i,] = unlist(temp.AVDS[1, ])
  record$sd.AVDS[i,] = unlist(temp.AVDS[2, ])
  
  if(i %*% 10 == 0){print(i)}
}


# Plot the results
result = lapply(record, FUN = function(x){apply(x, 2, mean)})

matplot(lambda.seq[-length(lambda.seq)], as.data.frame(result[c("oracle.infData", "oracle", "semiEfficient", "AVDS")]), type = "l", lwd = 2, lty = 1,
        # , ylim = c(-20, 20)
)
lines(lambda.seq[-length(lambda.seq)], result$semiEfficient + result$sd.semiEfficient, lty = 3, col = 3)
lines(lambda.seq[-length(lambda.seq)], result$semiEfficient - result$sd.semiEfficient, lty = 3, col = 3)
lines(lambda.seq[-length(lambda.seq)], result$AVDS + result$sd.AVDS, lty = 3, col = 5)
lines(lambda.seq[-length(lambda.seq)], result$AVDS - result$sd.AVDS, lty = 3, col = 5)
legend("topleft", legend = c("oracle.infData", "oracle", "semiEfficient", "AVDS"), lwd = 3, lty = 1, col = seq(1, 10))

# results
result
