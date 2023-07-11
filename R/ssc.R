#' Train an ssc model
#' @description Builds and trains a semi-supervised classifier with awareness of unknown classes, as described in Schrunner et al. (2020).
#' @param formula a formula object explaining the functional relation between input variables and target
#' @param data a data.frame containing the dataset
#' @export
ssc <- function(formula, data){
  res <- list()
  # implement training function

  class(res) <- "ssc"
  return(res)
}

print.ssc <- function(x){}

summary.ssc <- function(x){}

coef.ssc <- function(x){}

plot.ssc <- function(x){}

formula.ssc <- function(x){}

predict.ssc <- function(x){}

residuals.ssc <- function(x){}
