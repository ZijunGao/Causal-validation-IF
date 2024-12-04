# ihdp data
rm(list = ls())
source("~/Desktop/Research/Zijun/causal validation/Causal-validation/helper.R")

# install.packages("bartcs")
library("bartcs")
data = ihdp # 747 rows, 30 columns
data$y.0 = data$y_factual * (data$treatment == 0) + data$y_cfactual * (data$treatment == 1)
data$y.1 = data$y_cfactual * (data$treatment == 0) + data$y_factual * (data$treatment == 1)

# Preprocess the data
X.total = ihdp[, -seq(1, 5)]
n.total = dim(X.total)[1]
n = 300; n.train = n.total - n
d = dim(X.total)[2]

# Hyperparameter
n.fold = 2
method = "gradient boosting S learner"  # "LASSO", "gradient boosting S learner" 
prop.method = "LASSO" # "LASSO". The experiment is actually randomized.
sigma = 1 # Here we set sigma = 1 since sd((data$y_factual - data$mu1)[data$treatment == 1]), sd((data$y_factual - data$mu1)[data$treatment == 1]) ~ 1.

setting = method # "LASSO", "gradient boosting S learner" 

record.absolute.total = record.relative.total = list()
m = 100 # 100
record.absolute = record.relative = list()
record.absolute$oracle = record.absolute$semiEfficient = record.absolute$sd.oracle = record.absolute$sd.semiEfficient = record.absolute$AVDS = record.absolute$sd.AVDS = record.absolute$plug.in = record.absolute$sd.plug.in = matrix(0 , nrow = m, ncol = 2)
record.relative$oracle = record.relative$semiEfficient = record.relative$sd.oracle = record.relative$sd.semiEfficient = record.relative$AVDS = record.relative$sd.AVDS = matrix(0 , nrow = m, ncol = 1)

set.seed(200)
for(j in 1 : 1){ # only train once, estimate validation error m times.
  index = sample(n.total, n, replace = F)
  
  # training data
  X.train = as.matrix(X.total[-index,])  # covariates
  W.train = data$treatment[-index]  # binary treatment
  mu0.train = data$mu0[-index]
  mu1.train = data$mu1[-index]
  tau.train = mu1.train - mu0.train
  Y0.train = data$y.0[-index]
  Y1.train = data$y.1[-index]
  Y.train = Y1.train * W.train + Y0.train * (1 - W.train)
  
  # HTE estimator
  # LASSO
  lasso.fit.cv = cv.glmnet(x = cbind(X.train, W.train, W.train * X.train), Y.train, alpha = 1)
  
  # xgboost
  xgboost.fit = nuisance.learner(Y = Y.train, X = X.train, prop = mean(W.train), W = W.train, method = "gradient boosting S learner", train.index = seq(1, n.train), test.index = seq(1, n.train))
  
  # Define tau.hat.function based on the LASSO fit
  tau.hat.functions = function(X, j) {
    data.cnt = cbind(X, 0, 0 * X)
    data.trt = cbind(X, 1, 1 * X)
    if(j == "LASSO"){
      tau.hat = predict(lasso.fit.cv, newx = data.trt) - predict(lasso.fit.cv, newx = data.cnt)
    }else if(j == "xgboost"){
      data.new = data.frame(X, rep(0, dim(X)[1]), rep(0, dim(X)[1]))
      colnames(data.new) = c(colnames(X.train), "W", "Y")
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
  
  for(i in 1:m){
    # Validation data
    data$z = data$treatment
    data$y.0 = data$mu0 + rnorm(n.total, 0, sigma)
    data$y.1 = data$mu1 + rnorm(n.total, 0, sigma)
    
    X = as.matrix(X.total[index,])  # covariates
    W = data$z[index]  # binary treatment
    mu0 = data$mu0[index]
    mu1 = data$mu1[index]
    tau = mu1 - mu0
    Y0 = data$y.0[index]
    Y1 = data$y.1[index]
    Y = Y1 * W + Y0 * (1 - W)
    
    # evaluation
    folds = sample(rep(1 : n.fold, length.out = n))
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W,
                                 folds = folds, method = method, prop.method = prop.method) # debug linear
    
    # Absolute error
    # Oracle
    record.absolute$oracle[i, ] = sapply(c("LASSO", "xgboost"), function(x){mean((tau.hat.functions(X, x) - tau)^2)})
    record.absolute$sd.oracle[i, ] = sapply(c("xgboost"), function(x){sd((tau.hat.functions(X, x) - tau)^2 - (tau.hat.functions(X, "LASSO") - tau)^2) / sqrt(n)})  
    
    # Semi-parametrically efficient
    temp.semiEfficient = sapply(c("LASSO", "xgboost"), function(x){efficient.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
    record.absolute$semiEfficient[i,] = unlist(temp.semiEfficient[1,])
    record.absolute$sd.semiEfficient[i,] = unlist(temp.semiEfficient[2, ])
    
    # AVDS
    temp.AVDS = sapply(c("LASSO", "xgboost"), function(x){AVDS.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
    record.absolute$AVDS[i,] = unlist(temp.AVDS[1, ])
    record.absolute$sd.AVDS[i,] = unlist(temp.AVDS[2, ])
    
    # Plug-in
    temp.plug.in = sapply(c("LASSO", "xgboost"), function(x){plug.in.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, x)}, folds = folds, nuisance.hat = nuisance.hat)})
    record.absolute$plug.in[i,] = unlist(temp.plug.in[1, ])
    record.absolute$sd.plug.in[i,] = unlist(temp.plug.in[2, ])
    
    # Relative error
    record.relative$oracle[i] = sapply(c("xgboost"), function(x){mean((tau.hat.functions(X, x) - tau)^2)}) - mean((tau.hat.functions(X, "LASSO") - tau)^2)
    record.relative$sd.oracle[i] = sd((tau.hat.functions(X, "xgboost") - tau)^2- (tau.hat.functions(X, "LASSO") - tau)^2) / sqrt(n)
    
    # Semi-parametrically efficient
    temp.semiEfficient = efficient.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, "xgboost")}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, "LASSO")})
    record.relative$semiEfficient[i,] = unlist(temp.semiEfficient)[1]
    record.relative$sd.semiEfficient[i,] = unlist(temp.semiEfficient)[2]
    
    # AVDS
    temp.AVDS = AVDS.relative.MSE(Y = Y, X = X, W = W, tau.hat.function = function(X){tau.hat.functions(X, "xgboost")}, folds = folds, nuisance.hat = nuisance.hat, tau.hat.function.ref = function(X){tau.hat.functions(X, "LASSO")})
    record.relative$AVDS[i,] = unlist(temp.AVDS)[1]
    record.relative$sd.AVDS[i,] = unlist(temp.AVDS)[2]
  }
  record.absolute.total[[j]] = record.absolute
  record.relative.total[[j]] = record.relative
  print(j)
}

# Results
result = lapply(record.absolute.total, function(y){lapply(y, FUN = function(x){apply(x, 2, mean)})})
as.data.frame(result)

result = lapply(record.relative.total, function(y){lapply(y, FUN = function(x){mean(x)})})
as.data.frame(result)

hist(unlist(lapply(record.relative.total, function(x){mean(abs(x$AVDS - x$oracle) <= qnorm(0.95) * x$sd.AVDS)})), main = "AVDS", xlab = "", xlim = c(0, 1)); abline(v = 0.9, col = "red")
hist(unlist(lapply(record.relative.total, function(x){mean(abs(x$semiEfficient - x$oracle) <= qnorm(0.95) * x$sd.semiEfficient)})), main = "semiEfficient", xlab = "", xlim = c(0, 1)); abline(v = 0.9, col = "red")

# saveRDS(list(record.absolute.total = record.absolute.total, record.relative.total = record.relative.total), file.path("~/Desktop/Research/Zijun/causal validation/Causal-validation/data", paste("ihdp ", setting, ".rds", sep= "")))
