# convert formula and data.frame object to matrix X and vector y
convertFormula <- function(formula, data){
  frame <- model.frame(formula, data, na.action = na.pass)
  y <- as.factor(model.response(frame))
  formula <- update(formula, NULL ~ .) # remove target variable to keep rows in X with NAs in y
  X <- model.matrix(formula, data)
  return(list(X = X, y = y))
}

#' Train a Gaussian Bayes Classifier
#' @inheritParams SSC
#' @importFrom mvtnorm dmvnorm
BayesClassifier <- function(formula, data){
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  param <- list()
  logLik <- 0
  for(c in levels(y)){
    inds <- which(y == c)
    submatrix <- X[inds,,drop = F]

    mu <- colMeans(submatrix)
    Sigma <- cov(submatrix)
    prior <- length(inds) / length(y)

    logLik <- logLik + dmvnorm(submatrix, mu, Sigma, log = TRUE)
    param[[c]] <- list(mu = mu,
                       Sigma = Sigma,
                       prior = prior)
  }
  names(param) <- levels(y)

  model <- list(param = param,
                formula = formula,
                n = nrow(data),
                logLik = logLik)
  class(model) <- "BayesClassifier"
  return(model)
}

#' @importFrom mvtnorm dmvnorm
predict.BayesClassifier <- function(object, newdata, type = "probs"){
  l <- convertFormula(object$formula, newdata)
  X <- l$X

  # likelihood under class p
  probs <- c()
  for(p in object$param){
    probs <- cbind(probs,
                   dmvnorm(X,
                           p$mu,
                           p$Sigma,
                           log = TRUE) +
                    p$prior)
  }
  colnames(probs) <- names(object$param)

  if(type == "probs"){
    return(probs)
  }
  else{
    return(apply(probs, 1, which.max))
  }
}
