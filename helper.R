# To be deleted
# B = 2 * W * (W - prop.hat) / prop.hat / (1 - prop.hat)
# A = W - prop.hat
# IF5 = (1 - B) * (mu1.hat - mu0.hat)^2 + B * Y * (mu1.hat - mu0.hat) - A * (mu1.hat - mu0.hat)^2
# value5 = mean(IF5)
# 
# IF6 = (2 - B) * (mu1.hat - mu0.hat)^2
# value6 = mean(IF6)


# helper for causal validation
library("xgboost", "caret")
# Cross-fitting
# Input:
  # nuisance.learner: "linear regression"
# Example:
  # m = 1000; record = rep(0, m)
  # n = 400
  # folds = sample(rep(1 : 2, length.out = n))
  # cross_fitting_results = cross_fitting(Y = Y, X = X, W = W, folds = folds)
  # record[i] = second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = mu1.hat, mu0.hat = mu0.hat, prop.hat = p + 0.1) - mean(tau.test^2)
  # cross_fitting_results[folds]
cross_fitting = function(Y, X, W, folds = NULL, n.fold = 2, prop = NULL, method = "linear") {
  # preprocess
  n = length(Y)
  data = data.frame(Y, X, W, X * W)
  
  if(is.null(folds)){
    # split data into n.fold folds
    folds = sample(rep(1 : n.fold, length.out = n))
  }else{
    # number of folds
    n.fold = max(folds)
  }
  
  result = list()
  result$mu1.hat = result$mu0.hat = result$prop.hat = rep(0, n)
  for (k in 1 : n.fold) {
    # fit the model on the training data and predict on the test data here!!!
    if(is.null(prop)){prop.hat = mean(W[which(folds != k)])}else{prop.hat = prop}
    fit = nuisance.learner(Y = Y, X = X, prop = prop.hat, W = W, method = method, train.index = which(folds != k), test.index = which(folds == k))
    result$prop.hat[which(folds == k)] = rep(prop.hat, sum(folds == k))
    result$mu0.hat[which(folds == k)] = fit$mu0.hat
    result$mu1.hat[which(folds == k)] = fit$mu1.hat
  }
  return(result)
}

# nuisance learner
nuisance.learner = function(Y, X = NULL, prop = NULL, W = NULL, method = "linear", train.index = NULL, test.index = NULL, ...){
  n = length(Y)
  if(length(prop) == 1){prop = rep(prop, n)}
  
  if(method == "linear"){
    data.train = data.frame(Y, X, W - prop, (W-prop) * X)
    data.0 = data.frame(Y, X, 0 - prop, (0 - prop) * X); colnames(data.0) = colnames(data.train)
    data.1 = data.frame(Y, X, 1 - prop, (1 - prop) * X); colnames(data.1) = colnames(data.train)
    nuisance.model = lm(Y ~ ., data = data.train[train.index,])
    mu0.hat = predict(nuisance.model, newdata = data.0[test.index, ])
    mu1.hat = predict(nuisance.model, newdata = data.1[test.index, ])
    mu.hat = mu0.hat * (1 - prop[test.index]) + mu1.hat * prop[test.index]
    tau.hat = mu1.hat -  mu0.hat
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat)) 
    
  }else if(method == "gradient boosting"){
    data.full = data.frame(X, W - prop, Y)
    train.index_1 = createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data = data.full[train.index,]
    test.data = data.full[test.index,]
    
    num_cols = ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    # Fit the model with early stopping
    nrounds = 100 
    
    nuisance.mu = xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    # Make predictions on the training and validation data
    train.pred = predict(nuisance.mu, newdata = dtrain)
    
    # Calculate residuals for training and validation data
    train.residual = (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero
    
    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))
    
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)
    
    weights = abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error = function(preds, dtrain) {
      labels = getinfo(dtrain, "label")
      weights = getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad = weights * (preds - labels)
      hess = weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params = list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau = xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }
  else if(method == "gradient boosting early stopping"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 = createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data = data.full[train.index,][train.index_1, ]
    val.data = data.full[train.index,][-train.index_1, ]
    test.data = data.full[test.index,]
    
    num_cols = ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    #Watchlist to track performance on validation set
    watchlist = list(train = dtrain, eval = dval)
    
    # Fit the model with early stopping
    nrounds = 100  # Maximum number of boosting rounds
    early_stopping_rounds =10  # Stop early if there is no improvement
    
    
    nuisance.mu = xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    
    # Make predictions on the training and validation data
    train.pred = predict(nuisance.mu, newdata = dtrain)
    val.pred = predict(nuisance.mu, newdata = dval)

    # Calculate residuals for training and validation data
    train.residual = (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero
    val.residual = (val.data$Y - val.pred) / val.data[, num_cols-1]
    
    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))

    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.residual)
    
    watchlist = list(train = dtrain, eval = dval)
    weights = abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error = function(preds, dtrain) {
      labels = getinfo(dtrain, "label")
      weights = getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad = weights * (preds - labels)
      hess = weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params = list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau = xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }else if(method == "gradient boosting linear"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 = createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data = data.full[train.index,][train.index_1, ]
    val.data = data.full[train.index,][-train.index_1, ]
    test.data = data.full[test.index,]
    
    num_cols = ncol(data.full)
    features = 1:(num_cols-2)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    #Watchlist to track performance on validation set
    watchlist = list(train = dtrain, eval = dval)
    
    # Fit the model with early stopping
    nrounds = 100  # Maximum number of boosting rounds
    
    
    nuisance.mu = xgb.train(
      booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      # early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    
    # Make predictions on the training and validation data
    train.pred = predict(nuisance.mu, newdata = dtrain)
    val.pred = predict(nuisance.mu, newdata = dval)

    # Calculate residuals for training and validation data
    train.residual = (train.data$Y - train.pred) / train.data[, num_cols-1] # train.data[, num_cols-1] could be zero
    val.residual = (val.data$Y - val.pred) / val.data[, num_cols-1]
    
    
    dtest = xgb.DMatrix(data = as.matrix(test.data[,features]))
    
    
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.residual)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.residual)
    
    watchlist = list(train = dtrain, eval = dval)
    weights = abs(train.data[, num_cols-1])**2
    setinfo(dtrain, "weight", weights)
    
    weighted_squared_error = function(preds, dtrain) {
      labels = getinfo(dtrain, "label")
      weights = getinfo(dtrain, "weight")
      
      # Calculate the gradient and hessian
      grad = weights * (preds - labels)
      hess = weights
      
      return(list(grad = grad, hess = hess))
    }
    
    # Define the parameters for the second model
    params = list(eval_metric = "rmse")
    
    # Fit the second model with early stopping
    nuisance.tau = xgb.train(
      booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      # early_stopping_rounds = early_stopping_rounds,
      verbose = 0,
      obj = weighted_squared_error
    )
    
    #dtest.0 = xgb.DMatrix(data = as.matrix(data.0[test.index, -1]))
    #dtest.1 = xgb.DMatrix(data = as.matrix(data.1[test.index, -1]))
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat # 
    
    #mu0.hat = predict(nuisance.model, newdata = dtest.0)
    #mu1.hat = predict(nuisance.model, newdata = dtest.1)
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, mu.hat = mu.hat, tau.hat = tau.hat,tau=nuisance.tau)) 
  }
}

# Semi-parametrically efficient estimator of weighted ATE
# Input:
# Example:
  # weighted.ATE(Y = Y, X = X, W = W, mu1.hat = mu1, mu0.hat = mu0, prop.hat = p)
  # tau.hat = tau
  # mean(tau.hat^2) - 2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = mu1, mu0.hat = mu0, prop.hat = p, weight = tau) + mean(tau^2) 
weighted.ATE = function(Y, X, W, mu1.hat = 0, mu0.hat = 0, prop.hat = NULL, weight = NULL, IF = F){
  # preprocessf
  n = length(Y)
  if(is.null(weight)){weight = rep(1, n)}
  if(is.null(prop.hat)){prop.hat = rep(mean(W), n)}
  mu1.hat = c(mu1.hat)
  mu0.hat = c(mu0.hat)
  prop.hat = c(prop.hat)
  
  IF1 = weight * (W * (Y - mu1.hat) / prop.hat + mu1.hat - (1 - W) * (Y - mu0.hat) / (1 - prop.hat) - mu0.hat)
  value1 = mean(IF1)
  
  if(IF){return(IF1)}
  return(value1)
}

# Example
  # semi.parametric.relative.sd(Y = Y, X = X, W = W, mu1.hat = mu1, mu0.hat = mu0, tau.hat = tau.hat.function(X), tau.hat.alternative = tau.hat.function.alternative(X), prop.hat = p)
semi.parametric.relative.sd = function(Y, X, W,  mu1.hat,  mu0.hat, tau.hat = 0, tau.hat.alternative = 0, prop.hat = NULL){
  IF1 = tau.hat^2 - tau.hat.alternative^2 - 2 * (tau.hat - tau.hat.alternative) * (W * (Y - mu1.hat) / prop.hat + mu1.hat - (1 - W) * (Y - mu0.hat) / (1 - prop.hat) - mu0.hat)
  return(sd(IF1) / sqrt(length(Y)))
}

# semi-parametrically efficient estimator of second moment of ITE
# E[(Y(1) - Y(0) - tau.hat(X))^2] is not estimable, since it depends on both Y(1) and Y(0)
# Bias case: value 3 better
  # n = 400; d = 5; sigma = 1; p = 0.5
  # m = 1000; record = matrix(0, m, 6); beta0 = 0; beta = rep(1, d); delta = rep(1, d)
  # 
  # X.test = matrix(rnorm(n * 100 * d), nrow = n * 100)  # covariates
  # tau.test = (X.test %*% delta) + 2
  # for(i in 1:m){
  #   X = matrix(rnorm(n * d), nrow = n)  # covariates
  #   W = rbinom(n, 1, 0.5)  # binary treatment
  #   mu0 = beta0 + X %*% beta
  #   tau = (X %*% delta) + 2
  #   mu1 = mu0 + tau
  #   Y0 = mu0 + rnorm(n, 0, sigma)  # outcome with treatment effect
  #   Y1 = Y0 + tau
  #   Y = Y1 * W + Y0 * (1 - W)
  # 
  #   mu0.hat = beta0 + X %*% (beta + rnorm(d, 0, 5 / sqrt(n))) # mu0 + rnorm(n, 0, 5 / sqrt(n)) # unbiased is quite important
  #   mu1.hat = mu0.hat + X %*% (delta * 1 + rnorm(d, 0, 5 / sqrt(n))) + 2.2 # mu1 + rnorm(n, 0, 5 / sqrt(n)) + 2.2 # tau.hat should be incorrect
  #   record[i,] = second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = mu1.hat, mu0.hat = mu0.hat, prop.hat = p + 0.1) - mean(tau.test^2) # prop.hat = p + 0.1
  # }
  # data.frame(bias = apply(record, 2, mean),
  #            sd = apply(record, 2, sd))
second.moment.ITE = function(Y, X, W, mu1.hat = 0, mu0.hat = 0, prop.hat = 0.5, IF = F){
  # preprocess
  n = length(Y)
  if(is.null(prop.hat)){prop.hat = rep(mean(W), n)}
  mu1.hat = c(mu1.hat)
  mu0.hat = c(mu0.hat)
  prop.hat = c(prop.hat)

  IF =  (mu1.hat - mu0.hat)^2 + 2 * (mu1.hat - mu0.hat) * (W * (Y - mu1.hat) / prop.hat - (1 - W) * (Y - mu0.hat) / (1 - prop.hat))
  value3 = mean(IF)
  
  if(IF){return(IF3)}else{
    return(value3)
    # return(c(value1, value2, value3, value4, value5, value6))
  }
}

# confidence interval

absolute.error = function(X, Y, W, folds, n.fold = NULL, nuisance.hat, tau.test){
  # here!!!
}
