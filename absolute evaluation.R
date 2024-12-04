# causal validation
# absolute evaluation
  # use method = "gradient boosting S learner" and set n.rounds = 1
# TODO:
rm(list = ls())
source("~/Desktop/Research/Zijun/causal validation/Causal-validation/helper.R")

n = 2000 # number of observations; 400
n.train = 2000 # 400
n.test = 10000
d = 10  # number of covariates
n.fold = 2 # 2, 5

beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(1/4, d); # delta[1:(d/2)] = 0
sigma = 1 # 1
p = 0.5

method = "linear" # "linear"; "ridge"; "LASSO"; "gradient boosting"; "gradient boosting S learner"

m = 20
record = list()
record$oracle.infData = record$oracle = record$semiEfficient = record$sd.oracle.infData = record$sd.oracle = record$sd.semiEfficient = record$AVDS = record$sd.AVDS = record$plug.in = record$sd.plug.in = matrix(0 , nrow = m, ncol = 2)

setting = "confidence interval coverage" # "inaccurate nuisance function estimator"; "confidence interval coverage"

if(setting == "inaccurate nuisance function estimator"){
  n = 1000
  n.train = 1000
  m = 100 
  method = "gradient boosting S learner" 
}

# Test data
X.test = matrix(rnorm(n.test * d), nrow = n.test)  # covariates
W.test = rbinom(n.test, 1, p) # binary treatment
mu0.test = beta0 + X.test %*% beta # + X.test[, 1]^2
tau.test = X.test %*% delta + delta0 # + X.test[, 1]^2 
mu1.test = mu0.test + tau.test
Y0.test = mu0.test + rnorm(n.test, 0, sigma) # outcome with treatment effect
Y1.test = Y0.test + tau.test
Y.test = Y1.test * W.test + Y0.test * (1 - W.test)
for (i in 1 : m){
  # Generate data
  # training
  X.train = matrix(rnorm(n.train * d), nrow = n.train)  # covariates
  W.train = rbinom(n.train, 1, p) # binary treatment
  mu0.train = beta0 + X.train %*% beta # + X.train[, 1]^2
  tau.train = X.train %*% delta + delta0 # + X.train[, 1]^2 
  mu1.train = mu0.train + tau.train
  Y0.train = mu0.train + rnorm(n.train, 0, sigma)  # outcome with treatment effect
  Y1.train = Y0.train + tau.train
  Y.train = Y1.train * W.train + Y0.train * (1 - W.train)
  
  # validation
  X = matrix(rnorm(n * d), nrow = n)  # covariates
  W = rbinom(n, 1, 0.5)  # binary treatment
  mu0 = beta0 + X %*% beta # + X[, 1]^2
  tau = X %*% delta + delta0 # + X[, 1]^2
  mu1 = mu0 + tau
  Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  Y1 = Y0 + tau
  Y = Y1 * W + Y0 * (1 - W)
  
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
  
  # MSE
  folds = sample(rep(1 : n.fold, length.out = n))
  nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  
  # Oracle
  record$oracle.infData[i, ] = sapply(c("LASSO", "xgboost"), function(x){mean((tau.hat.functions(X.test, x) - tau.test)^2)})
  record$sd.oracle.infData[i, ] = sapply(c("xgboost"), function(x){sd((tau.hat.functions(X.test, x) - tau.test)^2 - (tau.hat.functions(X.test, "LASSO") - tau.test)^2) / sqrt(n.test)})
  
  record$oracle[i, ] = sapply(c("LASSO", "xgboost"), function(x){mean((tau.hat.functions(X, x) - tau)^2)})
  record$sd.oracle[i, ] = sapply(c("xgboost"), function(x){sd((tau.hat.functions(X, x) - tau)^2 - (tau.hat.functions(X, "LASSO") - tau)^2) / sqrt(n)})  
  
  # Semi-parametrically efficient
  temp.semiEfficient = sapply(c("LASSO", "xgboost"), function(x){efficient.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
  record$semiEfficient[i,] = unlist(temp.semiEfficient[1,])
  record$sd.semiEfficient[i,] = unlist(temp.semiEfficient[2, ])
  
  # AVDS
  temp.AVDS = sapply(c("LASSO", "xgboost"), function(x){AVDS.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
  record$AVDS[i,] = unlist(temp.AVDS[1, ])
  record$sd.AVDS[i,] = unlist(temp.AVDS[2, ])
  
  # Plug-in
  temp.plug.in = sapply(c("LASSO", "xgboost"), function(x){plug.in.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
  record$plug.in[i,] = unlist(temp.plug.in[1, ])
  record$sd.plug.in[i,] = unlist(temp.plug.in[2, ])
  
  if(i %*% 10 == 0){print(i)}
}

# results
result = lapply(record, FUN = function(x){apply(x, 2, mean)})
result

# confidence interval coverage
mean(abs(record$plug.in[, 1] - record$oracle[, 1]) <= qnorm(0.95) * record$sd.plug.in[, 1])
mean(abs(record$plug.in[, 2] - record$oracle[, 2]) <= qnorm(0.95) * record$sd.plug.in[, 2])
mean(abs(record$AVDS[, 1] - record$oracle[, 1]) <= qnorm(0.95) * record$sd.AVDS[, 1])
mean(abs(record$AVDS[, 2] - record$oracle[, 2]) <= qnorm(0.95) * record$sd.AVDS[, 2])
mean(abs(record$semiEfficient[, 1] - record$oracle[, 1]) <= qnorm(0.95) * record$sd.semiEfficient[, 1])
mean(abs(record$semiEfficient[, 2] - record$oracle[, 2]) <= qnorm(0.95) * record$sd.semiEfficient[, 2])

mean(abs(record$plug.in[, 1] - record$oracle.infData[, 1]) <= qnorm(0.95) * record$sd.plug.in[, 1])
mean(abs(record$plug.in[, 2] - record$oracle.infData[, 2]) <= qnorm(0.95) * record$sd.plug.in[, 2])
mean(abs(record$AVDS[, 1] - record$oracle.infData[, 1]) <= qnorm(0.95) * record$sd.AVDS[, 1])
mean(abs(record$AVDS[, 2] - record$oracle.infData[, 2]) <= qnorm(0.95) * record$sd.AVDS[, 2])
mean(abs(record$semiEfficient[, 1] - record$oracle.infData[, 1]) <= qnorm(0.95) * record$sd.semiEfficient[, 1])
mean(abs(record$semiEfficient[, 2] - record$oracle.infData[, 2]) <= qnorm(0.95) * record$sd.semiEfficient[, 2])

# pick the better HTE estimator
selection = list()
selection$oracle = 2 * (record$oracle[, 1] > record$oracle[, 2]) - 1
selection$oracle.infData = 2 * (record$oracle.infData[, 1] > record$oracle.infData[, 2]) - 1
selection$AVDS = ((record$AVDS[,1] - qnorm(0.975) * record$sd.AVDS[,1]) > (record$AVDS[,2] + qnorm(0.975) * record$sd.AVDS[,2])) -  ((record$AVDS[,2] - qnorm(0.975) * record$sd.AVDS[,2]) > (record$AVDS[,1] + qnorm(0.975) * record$sd.AVDS[,1])) 
selection$semiEfficient = ((record$semiEfficient[,1] - qnorm(0.975) * record$sd.semiEfficient[,1]) > (record$semiEfficient[,2] + qnorm(0.975) * record$sd.semiEfficient[,2])) -  ((record$semiEfficient[,2] - qnorm(0.975) * record$sd.semiEfficient[,2]) > (record$semiEfficient[,1] + qnorm(0.975) * record$sd.semiEfficient[,1])) 
lapply(selection, table)

table(selection$oracle, selection$AVDS)
table(selection$oracle, selection$semiEfficient)

# save results
# saveRDS(record, file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation/data", paste(setting, " absolute", ".rds", sep= "")))
