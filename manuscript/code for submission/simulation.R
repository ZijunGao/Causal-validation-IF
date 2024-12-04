# Simulated data
rm(list = ls())
source("helper.R")


n = 2000 # number of observations; 400
n.train = 2000 # 400
n.test = 10000
d = 10  # number of covariates
n.fold = 2 # 2, 5
beta0 = 1; beta = rep(1, d)
delta0 = 1; delta = rep(1, d)
sigma = 1 # 1
p = 0.5

method = "linear" # "linear"; "gradient boosting"; "gradient boosting early stopping", "gradient boosting S learner"

setting = "sample size" # "sample size", "nuisance learner"
if(setting == "sample size"){
  n.seq = c(500, 1000, 1500, 2000, 2500, 3000)
}else if(setting == "nuisance learner"){
  method.seq = c("true", "linear", "ridge", "LASSO", "gradient boosting", "gradient boosting S learner")  
}

# Hyperparameter
n.fold = 2
method = "linear" # "linear", "ridge", "LASSO", "gradient boosting S learner" 
prop.method = "linear" # "linear", "ridge", "LASSO", "gradient boosting"
sigma = 1

record.absolute.total.seq = record.relative.total.seq = list()
record.absolute.total = record.relative.total = list()
m = 100
record.absolute = record.relative = list()
record.absolute$oracle = record.absolute$semiEfficient = record.absolute$sd.oracle = record.absolute$sd.semiEfficient = record.absolute$AVDS = record.absolute$sd.AVDS = record.absolute$plug.in = record.absolute$sd.plug.in = matrix(0 , nrow = m, ncol = 2)
record.relative$oracle = record.relative$semiEfficient = record.relative$sd.oracle = record.relative$sd.semiEfficient = record.relative$AVDS = record.relative$sd.AVDS = matrix(0 , nrow = m, ncol = 1)

set.seed(200)
for(k in 1 : 6){
  if(setting == "sample size"){n = n.seq[k]
  }else if(setting == "nuisance learner"){method = method.seq[k]}
  for(j in 1 : 100){
    # training data
    X.train = matrix(rnorm(n.train * d), nrow = n.train)  # covariates
    colnames(X.train) = paste("x", seq(1, d), sep = "_")
    W.train = rbinom(n.train, 1, p) # binary treatment
    mu0.train = beta0 + X.train %*% beta # + X.train[, 1]^2
    tau.train = X.train %*% delta + delta0 # + X.train[, 1]^2 
    mu1.train = mu0.train + tau.train
    Y0.train = mu0.train + rnorm(n.train, 0, sigma)  # outcome with treatment effect
    Y1.train = Y0.train + tau.train
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
    
    for(i in 1 : m){
      # Validation data
      X = matrix(rnorm(n * d), nrow = n)  # covariates
      W = rbinom(n, 1, 0.5)  # binary treatment
      mu0 = beta0 + X %*% beta # + X[, 1]^2
      tau = X %*% delta + delta0 # + X[, 1]^2
      mu1 = mu0 + tau
      Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
      Y1 = Y0 + tau
      Y = Y1 * W + Y0 * (1 - W)
      
      # Evaluation
      folds = sample(rep(1 : n.fold, length.out = n))
      nuisance.hat = cross_fitting(Y = Y, X = X, W = W,
                                   folds = folds, method = method, prop.method = prop.method, mu0 = mu0, mu1 = mu1, prop = rep(p, n))
      
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
    print(c(k, j))
  }
  record.absolute.total.seq[[k]] = record.absolute.total
  record.relative.total.seq[[k]] = record.relative.total
}
# result = lapply(record.absolute, FUN = function(x){apply(x, 2, mean)})
# result 
# 
# data.frame(mean = unlist(lapply(record.relative, mean)),
#            sd = unlist(lapply(record.relative, sd)))
# 
# mean(abs(record.relative$AVDS - record.relative$oracle) <= qnorm(0.95) * record.relative$sd.AVDS)
# mean(abs(record.relative$semiEfficient - record.relative$oracle) <= qnorm(0.95) * record.relative$sd.semiEfficient)


result = lapply(record.absolute.total, function(y){lapply(y, FUN = function(x){apply(x, 2, mean)})})
as.data.frame(result)

# data.frame(mean = unlist(lapply(record.relative, mean)),
#            sd = unlist(lapply(record.relative, sd)))
result = lapply(record.relative.total, function(y){lapply(y, FUN = function(x){mean(x)})})
as.data.frame(result)

hist(unlist(lapply(record.relative.total, function(x){mean(abs(x$AVDS - x$oracle) <= qnorm(0.95) * x$sd.AVDS)})), main = "AVDS", xlab = "", xlim = c(0, 1)); abline(v = 0.9, col = "red")
hist(unlist(lapply(record.relative.total, function(x){mean(abs(x$semiEfficient - x$oracle) <= qnorm(0.95) * x$sd.semiEfficient)})), main = "semiEfficient", xlab = "", xlim = c(0, 1)); abline(v = 0.9, col = "red")
