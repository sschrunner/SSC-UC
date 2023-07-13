#' @import stats
convertFormula <- function(formula, data){
  frame <- model.frame(formula, data, na.action = na.pass)
  y <- as.factor(model.response(frame))
  formula <- update(formula, NULL ~ .) # remove target variable to keep rows in X with NAs in y
  X <- model.matrix(formula, data)
  return(list(X = X, y = y))
}

#' Bayes Classifier or Naive Bayes Classifier
#' @describeIn BayesClassifier train a BayesClassifier object
#' @inheritParams SSC
#' @importFrom mclust dmvnorm
#' @import stats
#' @export
BayesClassifier <- function(formula, data, naive = FALSE, var_eps = 1e-2){
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  param <- list()
  logLik <- 0
  for(c in levels(y)){
    inds <- which(y == c)
    submatrix <- X[inds,,drop = F]

    mu <- colMeans(submatrix)
    if(naive){
      Sigma <- diag(apply(submatrix, 2, var))
    }
    else{
      Sigma <- cov(submatrix)
    }
    Sigma <- Sigma + diag(rep(var_eps, ncol(X)))
    prior <- log(length(inds) / length(y))

    logLik <- logLik + sum(dmvnorm(submatrix, mu, Sigma, log = TRUE))
    param[[c]] <- list(mu = mu,
                       Sigma = Sigma,
                       prior = prior)
  }
  names(param) <- levels(y)

  model <- list(param = param,
                prior = prior,
                formula = formula,
                n = nrow(data),
                logLik = logLik,
                naive = naive)
  class(model) <- "BayesClassifier"
  return(model)
}

#' @describeIn BayesClassifier predict new data using a BayesClassifier object
#' @param object a BayesClassifier object
#' @param newdata a data.frame containing new features to predict
#' @param type type of return value, either "probs", or "class"
#' @param ... currently unused
#' @return either a matrix containing log-probabilities for each class, or a vector of classes, see 'type'.
#' @importFrom mclust dmvnorm
#' @export
predict.BayesClassifier <- function(object, newdata, type = "probs", ...){
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
  else if(type == "class"){
    return(colnames(probs)[apply(probs, 1, which.max)])
  }
  else{
    stop("Unknown type")
  }
}
