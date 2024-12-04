# causal validation
# absolute evaluation
# TODO:
rm(list = ls())
source("helper.R")

n = 1000  # number of observations; 1000 
n.train = 1000
n.test = 10000
d = 5  # number of covariates
n.fold = 2 # 2, 5

beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(1, d)
sigma = 1
p = 0.5 # 0.5

setting = "confidence interval coverage" # "inaccurate nuisance function estimator"; "confidence interval coverage"

method = "gradient boosting S learner" # "linear"; "gradient boosting"; "gradient boosting early stopping", "gradient boosting S learner"
m = 100
record = list()
record$oracle.infData = record$oracle = record$semiEfficient = record$sd.oracle.infData = record$sd.oracle = record$sd.semiEfficient = record$AVDS = record$sd.AVDS = matrix(0 , nrow = m, ncol = 1)

if(setting == "inaccurate nuisance function estimator"){
  n = 1000
  n.train = 1000
  m = 100 
  method = "gradient boosting S learner" 
}

set.seed(318)

# Test data
X.test = matrix(rnorm(n.test * d), nrow = n.test)  # covariates
W.test = rbinom(n.test, 1, p) # binary treatment
mu0.test = beta0 + X.test %*% beta # + X.test[, 1]^2
tau.test = X.test %*% delta + delta0 + X.test[, 1]^2 
mu1.test = mu0.test + tau.test
Y0.test = mu0.test + rnorm(n.test, 0, sigma)  # outcome with treatment effect
Y1.test = Y0.test + tau.test
Y.test = Y1.test * W.test + Y0.test * (1 - W.test)

# training
X.train = matrix(rnorm(n.train * d), nrow = n.train)  # covariates
W.train = rbinom(n.train, 1, p) # binary treatment
mu0.train = beta0 + X.train %*% beta # + X.train[, 1]^2
tau.train = X.train %*% delta + delta0 + X.train[, 1]^2 
mu1.train = mu0.train + tau.train
Y0.train = mu0.train + rnorm(n.train, 0, sigma)  # outcome with treatment effect
Y1.train = Y0.train + tau.train
Y.train = Y1.train * W.train + Y0.train * (1 - W.train)


# Fit LASSO with the current lambda
lasso.fit.cv = cv.glmnet(x = cbind(X.train, W.train, W.train * X.train), Y.train, alpha = 1)
# Fit xgboost with the current lambda
xgboost.fit = nuisance.learner(Y = Y.train, X = X.train, prop = mean(W.train), W = W.train, method = "gradient boosting S learner", train.index = seq(1, n.train), test.index = seq(1, n.train))

# Define tau.hat.function based on the LASSO fit
tau.hat.functions = function(X, j) {
  data.cnt = cbind(X, 0, 0 * X)
  data.trt = cbind(X, 1, 1 * X)
  if(j == "LASSO"){
    tau.hat = predict(lasso.fit.cv, newx = data.trt) - predict(lasso.fit.cv, newx = data.cnt)
  }else if(j == "xgboost"){
    data.new = data.frame(X, rep(0, dim(X)[1]), rep(0, dim(X)[1]))
    colnames(data.new) = c(paste("X", seq(1, dim(X)[2]), sep = ""), "W", "Y")
    new.data0 = data.new; new.data0$W = 0
    new.data1 = data.new; new.data1$W = 1
    
    dnew0 = xgb.DMatrix(data = as.matrix(new.data0[,seq(1, dim(X)[2]+1)]))
    dnew1 = xgb.DMatrix(data = as.matrix(new.data1[,seq(1, dim(X)[2]+1)]))
    mu0.hat = predict(xgboost.fit$nuisance.model, newdata = dnew0)
    mu1.hat = predict(xgboost.fit$nuisance.model, newdata = dnew1)
    tau.hat = mu1.hat - mu0.hat
  }
  return(tau.hat)
}
for (i in 1 : m){
  # validation
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, 0.5)  # binary treatment
  mu0 = beta0 + X %*% beta # + X[,1]^2 
  tau = X %*% delta + delta0 + X[,1]^2 # + X[,1]^2 
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
  # MSE
  folds = sample(rep(1 : n.fold, length.out = n))
  nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  
  # Oracle
  record$oracle.infData[i, ] = sapply(c("xgboost"), function(x){mean((tau.hat.functions(X.test, x) - tau.test)^2)}) - mean((tau.hat.functions(X.test, "LASSO") - tau.test)^2)
  record$sd.oracle.infData[i, ] = sapply(c("xgboost"), function(x){sd((tau.hat.functions(X.test, x) - tau.test)^2 - (tau.hat.functions(X.test, "LASSO") - tau.test)^2) / sqrt(n.test)})
  
  record$oracle[i, ] = sapply(c("xgboost"), function(x){mean((tau.hat.functions(X, x) - tau)^2)}) - mean((tau.hat.functions(X, "LASSO") - tau)^2)
  record$sd.oracle[i, ] = sapply(c("xgboost"), function(x){sd((tau.hat.functions(X, x) - tau)^2 - (tau.hat.functions(X, "LASSO") - tau)^2) / sqrt(n)})  
  # Semi-parametrically efficient
  temp.semiEfficient = efficient.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, "xgboost")}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, "LASSO")})
  record$semiEfficient[i,] = unlist(temp.semiEfficient)[1]
  record$sd.semiEfficient[i,] = unlist(temp.semiEfficient)[2]
  # AVDS
  temp.AVDS = AVDS.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, "xgboost")}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, "LASSO")})
  record$AVDS[i,] = unlist(temp.AVDS)[1]
  record$sd.AVDS[i,] = unlist(temp.AVDS)[2]
}

# results
data.frame(mean = unlist(lapply(record, mean)),
           sd = unlist(lapply(record, sd)))

# confidence interval coverage
mean(abs(record$AVDS - record$oracle) <= qnorm(0.95) * record$sd.AVDS)
mean(abs(record$semiEfficient - record$oracle) <= qnorm(0.95) * record$sd.semiEfficient)
# 
mean(abs(record$AVDS - record$oracle.infData) <= qnorm(0.95) * record$sd.AVDS)
mean(abs(record$semiEfficient - record$oracle.infData) <= qnorm(0.95) * record$sd.semiEfficient)

# pick the better HTE estimators
selection = list()
selection$oracle = ((record$oracle - qnorm(0.95) * record$sd.oracle) > 0) - ((record$oracle + qnorm(0.95) * record$sd.oracle) < 0)
selection$oracle.infData = ((record$oracle.infData - qnorm(0.95) * record$sd.oracle.infData) > 0) - ((record$oracle.infData + qnorm(0.95) * record$sd.oracle.infData) < 0)
selection$AVDS = ((record$AVDS - qnorm(0.95) * record$sd.AVDS) > 0) - ((record$AVDS + qnorm(0.95) * record$sd.AVDS) < 0)
selection$semiEfficient = ((record$semiEfficient - qnorm(0.95) * record$sd.semiEfficient) > 0) - ((record$semiEfficient + qnorm(0.95) * record$sd.semiEfficient) < 0)
lapply(selection, table)

table(selection$oracle, selection$AVDS)
table(selection$oracle, selection$semiEfficient)
