# helper for causal validation
library("xgboost", "caret")

# General input:
  # Y: A vector representing the responses.
  # X: A matrix representing the covariates. Each row corresponds to a data point.
  # W: A vector representing the treatment assignment.
  # mu0.hat: A vector of estimated E[Y(0) | X] based on the validation dataset, evaluated at X.
  # mu1.hat: A vector of estimated E[Y(1) | X] based on the validation dataset, evaluated at X.
  # prop.hat: A vector of estimated propensity scores E[W | X], evaluated at X.
  # tau.hat.function: The estimator of CATE to evaluate.
  # tau.hat.function.ref: The reference estimator of CATE is used when calculating relative performance. The relative error is determined by subtracting the error of the reference estimator tau.hat.function.ref from that of the estimator being evaluated tau.hat.function.
  # link.function: The function that converts the conditional mean to the natural parameter scale. The default is the identity link function: function(x){return(x)}.
  # link.function.derivative: The derivative function of the link.function. The default is the constant one function: function(x){return(rep(1, length(x)))}.
  # return.IF: Return the vector of IF plus the estimator based on estimated nuisance functions if TRUE. Otherwise, return the estimator value based on the IF, i.e., the mean of the IF.


# Cross-fitting
# Input:
# Example:
  # m = 1000; record = rep(0, m)
  # n = 400
  # folds = sample(rep(1 : 2, length.out = n))
  # cross_fitting_results = cross_fitting(Y = Y, X = X, W = W, folds = folds)
cross_fitting = function(Y, X, W, folds = NULL, n.fold = 2, prop = NULL, method = "linear", prop.method = NULL) {
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
    # fit the model on the training data and predict on the test data
    if(is.null(prop)){
      if(is.null(prop.method)){
        prop.hat = mean(W[which(folds != k)])
      }else{
        prop.hat = prop.nuisance.learner(X = X, W = W, method = prop.method, train.index = which(folds != k), test.index = which(folds == k)) 
      }
    }else{prop.hat = prop}
    fit = nuisance.learner(Y = Y, X = X, prop = prop.hat, W = W, method = method, train.index = which(folds != k), test.index = which(folds == k))
    result$prop.hat[which(folds == k)] = prop.hat[pmin(which(folds == k), length(prop.hat))] # length(prop.hat) may be one
    result$mu0.hat[which(folds == k)] = fit$mu0.hat
    result$mu1.hat[which(folds == k)] = fit$mu1.hat
  }
  return(result)
}

# nuisance learner
# Input:
  # train.index: Indices of the data points used for learning the nuisance functions.
  # test.index: Indices of the data points used for calculate the nuisance function values.
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
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat)) 
    
  }else if(method == "ridge"){
    data.1 = cbind(X, 1, 1 * X)
    data.0 = cbind(X, 1, 1 * X)
    nuisance.model = cv.glmnet(x = cbind(X, W, W * X)[train.index,], Y[train.index], alpha = 0)
    mu0.hat = predict(nuisance.model, newx = data.0[test.index, ])
    mu1.hat = predict(nuisance.model, newx = data.1[test.index, ])
    mu.hat = mu0.hat * (1 - prop[test.index]) + mu1.hat * prop[test.index]
    tau.hat = mu1.hat - mu0.hat
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat)) 
    
    }else if(method == "LASSO"){
      data.1 = cbind(X, 1, 1 * X)
      data.0 = cbind(X, 1, 1 * X)
      nuisance.model = cv.glmnet(x = cbind(X, W, W * X)[train.index,], Y[train.index], alpha = 1)
      mu0.hat = predict(nuisance.model, newx = data.0[test.index, ])
      mu1.hat = predict(nuisance.model, newx = data.1[test.index, ])
      mu.hat = mu0.hat * (1 - prop[test.index]) + mu1.hat * prop[test.index]
      tau.hat = mu1.hat - mu0.hat
      
      return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat)) 
      
    }else if(method == "gradient boosting"){
    data.full = data.frame(X, W - prop, Y)
    train.index_1 = caret::createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
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
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, tau=nuisance.tau)) 
  }
  else if(method == "gradient boosting early stopping"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 = caret::createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
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
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, tau=nuisance.tau)) 
  }else if(method == "gradient boosting linear"){
    
    data.full = data.frame(X, W - prop,Y)
    train.index_1 = caret::createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
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
    
    # Define the data matrix for the test data
    dtest = xgb.DMatrix(data = as.matrix(test.data[, 1:(num_cols-2)]))
    
    # Make predictions for the test data using the first and second models
    mu.hat = predict(nuisance.mu, newdata = dtest)
    tau.hat = predict(nuisance.tau, newdata = dtest)
    
    # Calculate the final predictions for mu0.hat and mu1.hat
    mu0.hat = mu.hat - prop[test.index] * tau.hat
    mu1.hat = mu0.hat + tau.hat # 
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, tau=nuisance.tau)) 
  }else if(method == "gradient boosting S learner"){    
    data.full = data.frame(X, W, Y)
    train.index_1 = caret::createDataPartition(data.full$Y[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data = data.full[train.index,][train.index_1, ]
    val.data = data.full[train.index,][-train.index_1, ]
    test.data = data.full[test.index,]
    test.data0 = test.data; test.data0$W = 0
    test.data1 = test.data; test.data1$W = 1
    
    num_cols = ncol(data.full)
    features = 1:(num_cols-1)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$Y)
    dval = xgb.DMatrix(data = as.matrix(val.data[, features]), label = val.data$Y)
    # Define the parameters
    params = list(objective = "reg:squarederror", eval_metric = "rmse")
    
    #Watchlist to track performance on validation set
    watchlist = list(train = dtrain, eval = dval)
    
    # Fit the model with early stopping
    nrounds = 100  # Maximum number of boosting rounds; 100; underfit: 10
    
    nuisance.model = xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      #lambda = 1,
      watchlist = watchlist,
      # early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    dtest0 = xgb.DMatrix(data = as.matrix(test.data0[,features]))
    dtest1 = xgb.DMatrix(data = as.matrix(test.data1[,features]))
    
    # Calculate the predictions for mu0.hat and mu1.hat
    mu0.hat = predict(nuisance.model, newdata = dtest0)
    mu1.hat = predict(nuisance.model, newdata = dtest1)
    
    return(result = list(mu0.hat = mu0.hat, mu1.hat = mu1.hat, nuisance.model = nuisance.model)) 
  }
}

# learner for propensity score
prop.nuisance.learner = function(X = NULL, W = NULL, method = "linear", train.index = NULL, test.index = NULL, ...){
  n = length(W)
  if(method == "linear"){
    data = data.frame(X, W)
    nuisance.model = glm(W ~ ., data = data[train.index,], family = "binomial")
    prop.hat = predict(nuisance.model, newdata = data, type = "response")
    return(prop.hat) 
  }else if(method == "ridge"){
    nuisance.model = cv.glmnet(x = X[train.index,], W[train.index], alpha = 0, family = "binomial")
    prop.hat = predict(nuisance.model, newx = X, type = "response")
    return(prop.hat) 
    
  }else if(method == "LASSO"){
    nuisance.model = cv.glmnet(x = X[train.index,], W[train.index], alpha = 1, family = "binomial")
    prop.hat = predict(nuisance.model, newx = X, type = "response")
    return(prop.hat) 
    
  }else if(method == "gradient boosting"){
    data.full = data.frame(X, W)
    train.index_1 = caret::createDataPartition(data.full$W[train.index], p = 0.9, list = FALSE) # 90% training data
    train.data = data.full[train.index,]
    test.data = data.full
    
    num_cols = ncol(data.full)
    features = 1:(num_cols-1)
    # Define the data matrices
    dtrain = xgb.DMatrix(data = as.matrix(train.data[, features]), label = train.data$W)
    # Define the parameters
    params = list(objective = "binary:logistic", eval_metric = "error")
    
    # Fit the model with early stopping
    nrounds = 100 
    
    nuisance.prop = xgb.train(
      #booster = "gblinear",
      params = params,
      data = dtrain,
      nrounds = nrounds,
      verbose = 0
    )
    
    # Make predictions on the training and validation data
    train.pred = predict(nuisance.prop, newdata = dtrain)
    
    # Make predictions for the test data using the first and second models
    dtest = xgb.DMatrix(data = as.matrix(test.data[, features]), label = test.data$W)
    
    prop.hat = predict(nuisance.prop, newdata = dtest)
    return(prop.hat)
  }
}


# Semi-parametrically efficient estimator 

# Semi-parametrically efficient estimator of weighted ATE
# Input:
  # weight: A vector of weights used by the weighted ATE E[weight * tau(X)].
# Example:
  # weighted.ATE(Y = Y, X = X, W = W, mu1.hat = mu1.hat, mu0.hat = mu0.hat, prop.hat = p)
weighted.ATE = function(Y, X, W, mu1.hat = 0, mu0.hat = 0, prop.hat = NULL, weight = NULL, link.function = NULL, link.function.derivative = NULL, return.IF = F){
  # preprocess
  if(is.null(weight)){weight = rep(1, n)}
  if(is.null(prop.hat)){prop.hat = rep(mean(W), length(W))}
  if(is.null(link.function)){
    link.function = function(x){return(x)}
    link.function.derivative = function(x){return(rep(1, length(x)))}
  }
  
  # Compute the IF plus the estimator based on estimated nuisance functions 
  IF = weight * (link.function.derivative(mu1.hat) * W * (Y - mu1.hat) / prop.hat + link.function(mu1.hat) - link.function.derivative(mu0.hat) * (1 - W) * (Y - mu0.hat) / (1 - prop.hat) - link.function(mu0.hat))

  if(return.IF){return(IF)}else{return(mean(IF))}
}


# Semi-parametrically efficient estimator of the second moment of CATE
# TODO:
  # Debug non-identity link.functions.
# Input: 
# Example:
  # n = 400; n.test = n * 100; d = 5; sigma = 1; p = 0.5
  # beta0 = 0; beta = rep(1, d); delta = rep(1, d)
  # m = 1000; record = rep(0, m)
  # 
  # X.test = matrix(rnorm(n.test * d), nrow = n.test)  # covariates
  # tau.test = (X.test %*% delta) + 2
  # for(i in 1 : m){
  #   X = matrix(rnorm(n * d), nrow = n)  # covariates
  #   W = rbinom(n, 1, 0.5)  # treatment
  #   mu0 = beta0 + X %*% beta
  #   tau = (X %*% delta) + 2
  #   mu1 = mu0 + tau
  #   Y0 = mu0 + rnorm(n, 0, sigma)
  #   Y1 = Y0 + tau
  #   Y = Y1 * W + Y0 * (1 - W)
  # 
  #   # nuisance function estimators
  #   # artificial
  #   mu0.hat = beta0 + X %*% (beta + rnorm(d, 0, 5 / sqrt(n)))
  #   mu1.hat = mu0.hat + X %*% (delta * 1 + rnorm(d, 0, 5 / sqrt(n))) + 2.5
  #   prop.hat = p + 0.1
  #   
  #   record[i] = second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = mu1.hat, mu0.hat = mu0.hat, prop.hat = prop.hat) - mean(tau.test^2)
  # }
  # data.frame(bias = mean(record),
  #            sd = sd(record))
second.moment.ITE = function(Y, X, W, mu1.hat = 0, mu0.hat = 0, prop.hat = 0.5, link.function = NULL, link.function.derivative = NULL, return.IF = F){
  # preprocess
  if(is.null(prop.hat)){prop.hat = rep(msean(W), length(W))}
  if(is.null(link.function)){
    link.function = function(x){return(x)}
    link.function.derivative = function(x){return(rep(1, length(x)))}
  }

  # Compute the IF plus the estimator based on estimated nuisance functions 
  IF = (link.function(mu1.hat) - link.function(mu0.hat))^2 + 2 * (link.function(mu1.hat)  - link.function(mu0.hat)) * (link.function.derivative(mu1.hat) * W * (Y - mu1.hat) / prop.hat - link.function.derivative(mu0.hat) * (1 - W) * (Y - mu0.hat) / (1 - prop.hat))
  
  if(return.IF){return(IF)}else{return(mean(IF))}
}


# Estimators for the MSE
# Semi-parametric efficient estimator for the MSE
# TODO:
  # Debug non-identity link.functions.
# Input:
  # n.fold: Number of folds for cross-fitting.
  # folds: Vector of fold membership, taking values from one to n.fold.
# Example:
  # tau.hat.function = function(X){return(X %*% beta)}
  # efficient.MSE(Y = Y, X = X, W = W, tau.hat.function = tau.hat.function)
efficient.MSE = function(X, Y, W, folds = NULL, n.fold = NULL, nuisance.hat = NULL, method = "linear", prop.method = NULL, tau.hat.function, link.function = NULL, link.function.derivative = NULL){
  if(is.null(n.fold)){
    if(is.null(folds)){n.fold = 2}else{n.fold = max(folds)}
  }
  if(is.null(folds)){folds = sample(rep(1 : n.fold, length.out = n))}
  if(is.null(nuisance.hat)){
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method, prop.method = prop.method)
  }
  
  # Semi-efficient estimator based on the derived EIF
  record = list()
  IF.semiEfficient = tau.hat.function(X)^2 + second.moment.ITE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, prop.hat = nuisance.hat$prop.hat, link.function = link.function, link.function.derivative = link.function.derivative, return.IF = TRUE) - 
  2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, prop.hat = nuisance.hat$prop.hat, weight = tau.hat.function(X), link.function = link.function, link.function.derivative = link.function.derivative, return.IF = TRUE)
  record$error = mean(IF.semiEfficient)
  record$sd = sd(IF.semiEfficient) / sqrt(length(Y))
  
  return(record)
}


# Estimator for the MSE in the paper by Alaa et. al.
# TODO:
  # Debug non-identity link.functions.
# Input:
# Example:
AVDS.MSE = function(X, Y, W, folds = NULL, n.fold = NULL, nuisance.hat = NULL, method = "linear", tau.hat.function){
  if(is.null(n.fold)){
    if(is.null(folds)){n.fold = 2}else{n.fold = max(folds)}
  }
  if(is.null(folds)){folds = sample(rep(1 : n.fold, length.out = n))}
  if(is.null(nuisance.hat)){
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  }
  
  record = list()
  B = 2 * W * (W - nuisance.hat$prop.hat) / nuisance.hat$prop.hat / (1 - nuisance.hat$prop.hat)
  A = W - nuisance.hat$prop.hat
  IF = (1 - B) * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)^2 + 
    B * Y * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X)) - 
    A * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X))^2 + tau.hat.function(X)^2
  record$error = mean(IF)
  record$sd = sd(IF) / sqrt(n)
  
  return(record)
}


# Plug-in estimator for the MSE
# TODO:
  # Debug non-identity link.functions.
# Input:
# Example:
plug.in.MSE = function(X, Y, W, folds = NULL, n.fold = NULL, nuisance.hat = NULL, method = "linear", tau.hat.function, link.function = NULL, link.function.derivative = NULL){
  if(is.null(link.function)){
    link.function = function(x){return(x)}
    link.function.derivative = function(x){return(rep(1, length(x)))}
  }
  if(is.null(n.fold)){
    if(is.null(folds)){n.fold = 2}else{n.fold = max(folds)}
  }
  if(is.null(folds)){folds = sample(rep(1 : n.fold, length.out = n))}
  if(is.null(nuisance.hat)){
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  }
  
  record = list()
  IF = (tau.hat.function(X) - (link.function.derivative(nuisance.hat$mu1.hat) * W * (Y - nuisance.hat$mu1.hat) / nuisance.hat$prop.hat + link.function(nuisance.hat$mu1.hat) - link.function.derivative(nuisance.hat$mu0.hat) * (1 - W) * (Y - nuisance.hat$mu0.hat) / 
                            (1 - nuisance.hat$prop.hat) - link.function(nuisance.hat$mu0.hat)))^2
  record$error = mean(IF)
  record$sd = sd(IF) / sqrt(n)
  
  return(record)
}


# Semi-parametrically efficient estimator of the relative difference of the second moment of CATE
# Input:
# Example:
efficient.relative.MSE = function(X, Y, W, folds = NULL, n.fold = NULL, nuisance.hat = NULL, method = "linear", tau.hat.function, tau.hat.function.ref, link.function = NULL, link.function.derivative = NULL){
  if(is.null(n.fold)){
    if(is.null(folds)){n.fold = 2}else{n.fold = max(folds)}
  }
  if(is.null(folds)){folds = sample(rep(1 : n.fold, length.out = n))}
  if(is.null(nuisance.hat)){
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  }
  
  # Semi-efficient estimator based on the derived EIF
  record = list()
  IF.semiEfficient = tau.hat.function(X)^2 - tau.hat.function.ref(X)^2 -
    2 * weighted.ATE(Y = Y, X = X, W = W, mu1.hat = nuisance.hat$mu1.hat, mu0.hat = nuisance.hat$mu0.hat, prop.hat = nuisance.hat$prop.hat, weight = tau.hat.function(X) - tau.hat.function.ref(X), link.function = link.function, link.function.derivative = link.function.derivative, return.IF = TRUE)
  record$error = mean(IF.semiEfficient)
  record$sd = sd(IF.semiEfficient) / sqrt(length(Y))
  
  return(record)
}


# Estimator of the relative difference of the second moment of CATE based on the paper by Alaa et. al.
# Input:
# Example:
AVDS.relative.MSE = function(X, Y, W, folds = NULL, n.fold = NULL, nuisance.hat = NULL, method = "linear", tau.hat.function, tau.hat.function.ref){
  if(is.null(n.fold)){
    if(is.null(folds)){n.fold = 2}else{n.fold = max(folds)}
  }
  if(is.null(folds)){folds = sample(rep(1 : n.fold, length.out = n))}
  if(is.null(nuisance.hat)){
    nuisance.hat = cross_fitting(Y = Y, X = X, W = W, folds = folds, method = method)
  }
  
  record = list()
  B = 2 * W * (W - nuisance.hat$prop.hat) / nuisance.hat$prop.hat / (1 - nuisance.hat$prop.hat)
  A = W - nuisance.hat$prop.hat
  IF = (1 - B) * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)^2 + 
    B * Y * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X)) - 
    A * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function(X))^2 + tau.hat.function(X)^2
  IF.ref = (1 - B) * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat)^2 + 
    B * Y * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function.ref(X)) - 
    A * (nuisance.hat$mu1.hat - nuisance.hat$mu0.hat - tau.hat.function.ref(X))^2 + tau.hat.function.ref(X)^2
  record$error = mean(IF - IF.ref)
  record$sd = sd(IF - IF.ref) / sqrt(n)
  
  return(record)
}

