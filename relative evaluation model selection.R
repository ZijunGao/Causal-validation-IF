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
library(glmnet)

# Modular function to calculate errors given a specific tau.hat function and a fixed dataset
absolute_errors = function(tau.hat.function, X, W, Y, tau, nuisance.hat) {
  record = list()
  n = length(Y)
  
  # Oracle error
  record$error.oracle = mean((tau.hat.function(X) - tau)^2)
  
  # Semi-efficient estimator based on the derived EIF
  IF.semiEfficient = tau.hat.function(X)^2 + 
    second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, 
                      prop.hat = nuisance.hat$prop.hat, IF = TRUE) - 
    2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, 
                     prop.hat = nuisance.hat$prop.hat, weight = tau.hat.function(X), IF = TRUE)
  record$error.semiEfficient = mean(IF.semiEfficient)
  record$sd.semiEfficient = sd(IF.semiEfficient) / sqrt(n)
  
  # Theorem 2 method from the paper
  B = 2 * W * (W - nuisance.hat$prop.hat) / nuisance.hat$prop.hat / (1 - nuisance.hat$prop.hat)
  A = W - nuisance.hat$prop.hat
  IF.IF = (1 - B) * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)^2 + 
    B * Y * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X)) - 
    A * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X))^2 + tau.hat.function(X)^2
  record$error.IF = mean(IF.IF)
  record$sd.IF = sd(IF.IF) / sqrt(n)
  
  # Baseline method
  IF.baseline = (tau.hat.function(X) - (W * (Y - nuisance.hat$mu1.hat) / nuisance.hat$prop.hat + 
                                          nuisance.hat$mu1.hat - (1 - W) * (Y - nuisance.hat$mu0.hat) / 
                                          (1 - nuisance.hat$prop.hat) - nuisance.hat$mu0.hat))^2
  record$error.baseline = mean(IF.baseline)
  
  return(record)
}


# Main function to apply the sequence of tau.hat functions to the same dataset
absolute_errors_multiple_tau_hat = function(Y, X, W, tau, tau.hat.functions, n.fold = 5, method = "linear") {
  n = length(Y)
  d = dim(X)[2]
  
  # Estimate nuisance functions via cross-fitting
  folds = sample(rep(1:n.fold, length.out = n))
  nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  
  # Apply each tau.hat function to the same data
  results = list()
  
  # LASSO-based evaluation
  for (j in 1:length(lambda.seq)) {
    results[[j]] = absolute_errors(tau.hat.function = function(X.new) {return(tau.hat.functions(X.new)[, j])}, X, W, Y, tau, nuisance.hat)
  }
  
  results = do.call(rbind, lapply(results, as.data.frame))
  return(results)
}

# Apply the sequence to the same dataset and report errors
n = 400; n.test = 1000; d = 5; n.fold = 2 # 5; 2
beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(0, d); delta[1:2] = rep(2, 2)
sigma = 1
p = 0.5 # 0.5; 0.2

lambda.seq = 2^(-seq(-1, 4) / 1) # simple: 2^(-seq(-1, 4) / 1); complex: 2^(-seq(-5, 20) / 5)
m = 100
record = array(0, dim = c(m, length(lambda.seq), 6))

X.test = matrix(rnorm(n.test * d), nrow = n.test)  # covariates
W.test = rbinom(n.test, 1, p) # binary treatment
mu0.test = beta0 + X.test %*% beta # + X[, 1]^2
tau.test = X.test %*% delta + delta0 # + X[, 1]^2 
mu1.test = mu0.test + tau.test
Y0.test = mu0.test + rnorm(n.test, 0, sigma)  # outcome with treatment effect
Y1.test = Y0.test + tau.test
Y.test = Y1.test * W.test + Y0.test * (1 - W.test)

for(i in 1 : m){
  # Generate  data
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, p) # binary treatment
  mu0 = beta0 + X %*% beta # + X[, 1]^2
  tau = X %*% delta + delta0 # + X[, 1]^2 
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
  # Fit LASSO with the current lambda
  lasso.fit = glmnet(x = cbind(X, W, W * X), Y, alpha = 1, lambda = lambda.seq)
  
  # Define tau.hat.function based on the LASSO fit
  tau.hat.functions = function(X) {
    data.cnt = cbind(X, 0, 0 * X)
    data.trt = cbind(X, 1, 1 * X)
    tau.hat = predict(lasso.fit, newx = data.trt) - predict(lasso.fit, newx = data.cnt)
    return(tau.hat)
  }
  
  # Get results from the LASSO tau.hat functions
  record[i,,] = as.matrix(absolute_errors_multiple_tau_hat(Y = Y.test, X = X.test, W = W.test, tau = tau.test, n.fold = n.fold, tau.hat.functions = tau.hat.functions))
  
  if(i %*% 10 == 0){print(i)}
}

# Plot the results
matplot(lambda.seq, apply(record[, ,c(1, 2, 4, 6)], c(2, 3), mean), type = "l", lwd = 2, col = c(1, 2, 3, 4)
        # , ylim = c(-20, 20)
)
lines(lambda.seq, apply(record[, ,2], 2, mean) + apply(record[,,3], 2, mean), lty = 3, col = "red")
lines(lambda.seq, apply(record[, ,2], 2, mean) - apply(record[,,3], 2, mean), lty = 3, col = "red")
