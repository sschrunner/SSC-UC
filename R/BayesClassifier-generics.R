#' BayesClassifier generics
#' @description Generic functions implemented for the BayesClassifier class
#' @param x a BayesClassifier object
#' @param ... currently unused
#' @describeIn print.BayesClassifier prints the BayesClassifier
#' @returns The form of the returned values of the generic functions for class ´BayesClassifier´ are as follows:
#' \itemize{
#'  \item{´print´ and ´summary´ return a console output only}
#'  \item{´coef´ returns the list of parameters, see \link{BayesClassifier}}
#'  \item{´dimnames´ and ´levels´ return a vector containing the feature names and class names, respectively}
#'  \item{´dim´, ´length´ and ´nobs´ return a scalar representing the number of features, the number of classes and the number of training observations, respectively}
#'  \item{´formula´ returns the formula object used to call the BayesClassifier}
#'  \item{´logLik´, ´AIC´ and ´BIC´ return the log-likelihood, the AIC and the BIC values, respectively}
#' }
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
#' @importFrom utils head
#' @export
summary.BayesClassifier <- function(object, ...){
  cat("BayesClassifier model with", length(object), "classes and", dim(object), "non-constant features\n")
  if(dim(object) < length(object$all.features)){
    cat("Note: the total number of features is", length(object$all.features), "\n")
  }
  cat(rep("=",30), "\n", sep = "")
  cat("formula: ", deparse(object$formula), "\n")
  cat("used features: ", paste0(dimnames(object), collapse = ", "), "\n")
  cat("parameters: ")
  df <- data.frame(
       mu = sapply(coef(object),
                   function(x){return(paste0(round(head(x$mu, n = 10), 2), collapse = ","))}),
       Sigma = sapply(coef(object),
                      function(x){return(paste0(round(head(as.vector(x$Sigma), n = 10), 2), collapse = ","))}),
       prior = sapply(coef(object),
                      function(x){return(ifelse(exp(x$prior) > 0.005, round(exp(x$prior),2), "<0.01"))})
  )
  if(dim(object) > 10){
    df$mu <- paste0(df$mu, ",...")
    df$Sigma <- paste0(df$Sigma, ",...")
  }
  kable(df)
}

#' @describeIn print.BayesClassifier returns the list of model parameters in the BayesClassifier
#' @export
coef.BayesClassifier <- function(object, ...){
  return(object$param)
}

#' @describeIn print.BayesClassifier returns the names the features used in the BayesClassifier
#' @export
dimnames.BayesClassifier <- function(x){
  return(names(x$param[[1]]$mu))
}

#' @describeIn print.BayesClassifier returns the dimensionality of the feature space used in the BayesClassifier
#' @export
dim.BayesClassifier <- function(x){
  return(length(x$param[[1]]$mu))
}

#' @describeIn print.BayesClassifier returns the names of the classes in the BayesClassifier
#' @export
levels.BayesClassifier <- function(x){
  return(names(x$param))
}

#' @describeIn print.BayesClassifier returns the number of classes in the BayesClassifier
#' @export
length.BayesClassifier <- function(x){
  return(length(levels(x)))
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
