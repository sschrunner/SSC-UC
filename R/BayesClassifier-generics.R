#' BayesClassifier generics
#' @description Generic functions implemented for the BayesClassifier class
#' @param x a BayesClassifier object
#' @param ... currently unused
#' @describeIn print.BayesClassifier prints the BayesClassifier
#' @importFrom utils str
#' @export
print.BayesClassifier <- function(x, ...){
  print(x$formula)
  print(str(x$param))
}

#' @describeIn print.BayesClassifier summarizes the BayesClassifier
#' @param object a BayesClassifier object
#' @importFrom stats coef
#' @importFrom knitr kable
#' @export
summary.BayesClassifier <- function(object, ...){
  cat("BayesClassifier model with ", dim(object), " classes \n")
  print(object$formula)

  cat("Model parameters")
  kable(
    data.frame(
      mu = sapply(coef(object), function(x){return(paste0(round(x$mu, 2), collapse = ","))}),
      Sigma = sapply(coef(object), function(x){return(paste0(round(x$Sigma, 2), collapse = ","))}),
      prior = sapply(coef(object), function(x){return(ifelse(exp(x$prior) > 0.005, round(exp(x$prior),2), "<0.01"))})
    )
  )
}

#' @describeIn print.BayesClassifier returns the list of model parameters in the BayesClassifier
#' @export
coef.BayesClassifier <- function(object, ...){
  return(object$param)
}

#' @describeIn print.BayesClassifier returns the number of classes in the BayesClassifier
#' @export
dim.BayesClassifier <- function(x){
  return(length(x$param))
}

#' @describeIn print.BayesClassifier returns the number of training observations used for the BayesClassifier
#' @export
nobs.BayesClassifier <- function(object, ...){
  return(object$n)
}

#' @describeIn print.BayesClassifier returns the formula of the BayesClassifier
#' @export
formula.BayesClassifier <- function(x, ...){
  return(x$formula)
}

#' @describeIn print.BayesClassifier returns the log-likelihood of the BayesClassifier
#' @export
logLik.BayesClassifier <- function(object, ...){
  return(object$logLik)
}

#' @describeIn print.BayesClassifier returns the Akaike Information Criterion for the BayesClassifier
#' @param k penalty term; k=2 in classical AIC
#' @importFrom stats AIC
#' @export
AIC.BayesClassifier <- function(object, ..., k = 2){
  num_param <- ifelse(object$naive,
                      dim(object) * (2 * (-0.5 + sqrt(length(unlist(coef(object))) / dim(object) - 0.75)) + 1),
                      length(unlist(coef(object))))
  return(k * num_param - 2 * logLik(object))
}

#' @describeIn print.BayesClassifier returns the Bayesian Information Criterion for the BayesClassifier
#' @importFrom stats BIC
#' @export
BIC.BayesClassifier <- function(object, ...){
  return(AIC(object, k = log(nobs(object)), ...))
}
