#' @import stats
convertFormula <- function(formula, data){
  frame <- model.frame(formula, data, na.action = "na.pass")
  y <- as.factor(model.response(frame))
  X <- model.matrix(formula, frame)
  return(list(X = X, y = y))
}

#' Bayes Classifier or Naive Bayes Classifier
#' @description Train a BayesClassifier object
#' @details The function trains a Bayes classifier on a given classification dataset, specified by a data.frame ´data´ and a formula ´formula´. The specified target variable must be a factor.
#' The function estimates class-wise mean values (´mu´), covariance matrices (´Sigma´) and prior probabilities (´prior´) to represent classes.
#' If a Naive Bayes classifier is selected, only a diagonal covariance matrix is estimated.
#' Prior options indicate whether a uniform prior (all classes are equally weighted), or a proportional prior (all classes are weighted by their proportions in the training data) should be used.
#' @inheritParams SSC
#' @returns an S3 object of class ´BayesClassifier´, which has the following internal structure:
#' \itemize{
#'  \item{´param´: a list of model parameters; each list element represents one class, and contains a vector mu (class mean), a matrix Sigma (class covariance), and a scalar prior (prior probability)}
#'  \item{´prior´: the type of prior model used to call the BayesClassifier}
#'  \item{´formula´: the formula used to call the BayesClassifier}
#'  \item{´naive´: whether the model contains a Naive Bayes or a full Bayes classifier}
#'  \item{´n´: the number of samples during training}
#'  \item{´all.features´: the vector of all feature names in the training data}
#'  \item{´logLik´: the value of the log-likelihood}
#' }
#' @seealso \link{SSC} for semi-supervised classification with known and unknown classes, \link{EM} for semi-supervised classification with known classes only
#' @examples
#'  set.seed(1)
#'  dat <- data.frame(y = rep(c(1,2,3), each = 5),
#'     x1 = as.vector(sapply(c(5,0,0), rnorm, n = 5)),
#'     x2 = as.vector(sapply(c(0,0,5), rnorm, n = 5)))
#'  mod <- BayesClassifier(y ~ ., dat, naive = TRUE)
#'  summary(mod)
#'  predict(mod, newdata = expand.grid(x1 = c(0,5), x2 = c(0,5), y = NA), type = "class")
#' @importFrom mclust dmvnorm
#' @import stats
#' @export
BayesClassifier <- function(formula, data,
                            naive = FALSE, prior = "proportional",
                            var_eps = 1e-2){

  # convert formula to X,y format
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  # drop unused levels & remove data points without label
  y_nas <- is.na(y)
  if(any(y_nas)){
    warning(paste("BayesClassifier removed", sum(y_nas), "NAs"))
  }
  X <- X[!y_nas,]
  y <- y[!y_nas]
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
    prior_prob <- ifelse(prior == "proportional",
                         log(length(inds) / length(y)),
                         log(1 / nlevels(y)))

    # update log-likelihood
    logLik <- logLik + sum(dmvnorm(submatrix, mu, Sigma, log = TRUE))

    param[[c]] <- list(mu = mu,
                       Sigma = Sigma,
                       prior = prior_prob)
  }
  names(param) <- levels(y)

  model <- list(param = param,
                prior = prior,
                formula = formula,
                naive = naive,
                n = nrow(data),
                all.features = all_feats,
                logLik = logLik)
  class(model) <- "BayesClassifier"
  return(model)
}

#' Predictions from a BayesClassifier model
#' @description Predicts new data using a BayesClassifier object
#' @details Given a BayesClassifier object and a data.frame of new data, predictions are computed for each data point.
#' Two types of predictions are available:
#' \itemize{
#'  \item{´probs´: returns the matrix of log-class probabilities for each sample and class}
#'  \item{´class´: returns the vector of class labels with the highest predictive probabilities for each sample}
#' }
#' @param object a BayesClassifier object
#' @param newdata a data.frame containing new features to predict
#' @param type type of return value, either "probs", or "class"
#' @param ... currently unused
#' @return either a matrix containing log-probabilities, or a vector of classes, see ´details´.
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
