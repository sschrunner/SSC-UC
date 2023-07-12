#' @export
print.BayesClassifier <- function(x, ...){
  print(x$formula)
  print(str(x$param))
}

#' @importFrom knitr kable
#' @export
summary.BayesClassifier <- function(object, ...){
  cat("BayesClassifier model with ", dim(object), " classes \n")
  print(formula)

  cat("Model parameters")
  kable(
    data.frame(
      mu = sapply(coef(ssc), function(x){return(paste0(round(x$mu, 2), collapse = ","))}),
      Sigma = sapply(coef(ssc), function(x){return(paste0(round(x$Sigma, 2), collapse = ","))}),
      prior = sapply(coef(ssc), function(x){return(ifelse(x$prior > 0.005, round(x$prior,2), "<0.01"))})
    )
  )
}

#' @export
coef.BayesClassifier <- function(object, ...){
  return(object$param)
}

#' @export
dim.BayesClassifier <- function(x){
  return(length(x$param))
}

#' @export
formula.BayesClassifier <- function(x, ...){
  return(x$formula)
}

#plot.BayesClassifier <- function(x){}

#' @export
logLik.BayesClassifier <- function(object, ...){
  return(object$logLik)
}

#' @export
BIC.BayesClassifier <- function(object, ...){
  return(dim(object) * log(object$n) - 2 * logLik(object))
}
