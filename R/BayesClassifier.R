#' @import stats
convertFormula <- function(formula, data){
  frame <- model.frame(formula, data, na.action = "na.pass")
  y <- as.factor(model.response(frame))
  X <- model.matrix(formula, frame)
  return(list(X = X, y = y))
}

#' Bayes Classifier or Naive Bayes Classifier
#' @describeIn BayesClassifier train a BayesClassifier object
#' @inheritParams SSC
#' @examples
#'  set.seed(1)
#'  dat <- data.frame(y = rep(c(1,2,3), each = 5),
#'     x1 = as.vector(sapply(c(5,0,0), rnorm, n = 5)),
#'     x2 = as.vector(sapply(c(0,0,5), rnorm, n = 5)))
#'  mod <- BayesClassifier(y ~ x1 + x2, dat, naive = TRUE)
#'  summary(mod)
#'  predict(mod, newdata = expand.grid(x1 = c(0,5), x2 = c(0,5), y = NA), type = "class")
#' @importFrom mclust dmvnorm
#' @import stats
#' @export
BayesClassifier <- function(formula, data, naive = FALSE, var_eps = 1e-2){

  # convert formula to X,y format
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  # drop unused levels
  y <- droplevels(y)

  # determine columns with constant values across all classes
  all_feats <- colnames(X)
  const_cols <- apply(X, 2, function(x){return(length(unique(x)))}) == 1
  X <- X[,!const_cols]

  # initialize param and log-likelihood
  param <- list()
  logLik <- 0

  # estimate parameters for each level
  for(c in levels(y)){
    # extract submatrix of labeled data with level c
    inds <- which(y == c)
    submatrix <- X[inds,,drop = F]

    # estimate mu and Sigma for class c
    mu <- colMeans(submatrix)
    if(naive){
      # only variances for naive Bayes
      Sigma <- diag(apply(submatrix, 2, var))
    }
    else{
      # full covariance matrix for Bayes
      Sigma <- cov(submatrix)
    }
    # add constant variance term to avoid 0-entries
    Sigma <- Sigma + diag(rep(var_eps, ncol(X)))

    # estimate prior
    prior <- log(length(inds) / length(y))

    # update log-likelihood
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
                all.features = all_feats,
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

  # subset X to non-zero dimensions used by classifier
  non_const_cols <- dimnames(object)
  X <- X[,non_const_cols]

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
