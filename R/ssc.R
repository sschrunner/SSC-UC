#' Train an SSC-UC model
#' @description Builds and trains a semi-supervised classifier with awareness of unknown classes, as described in Schrunner et al. (2020).
#' @param formula a formula object explaining the functional relation between input variables and target
#' @param data a data.frame containing the dataset
#' @param perc_spies percentage of unlabeled points sampled as spies
#' @param naive whether a Bayes or Naive Bayes classifier should be used
#' @param var_eps scalar to add to the main diagonal of the covariance matrix to assure numerical stability. If 0, no scalar will be added.
#' @param max_unknown maximum number of unknown classes
#' @return a BayesClassifier object
#' @import mclust
#' @import stats
#' @export
SSC <- function(formula, data,
                perc_spies = 0.05, naive = FALSE, var_eps = 1e-2, max_unknown = 6){
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  # extract information on labeled and unlabeled samples
  labeled_inds <- which(!is.na(y))
  num_labeled <- length(labeled_inds)
  unlabeled_inds <- which(is.na(y))
  num_unlabeled <- length(unlabeled_inds)

  # implement training function
  S <- sample(labeled_inds, size = round(perc_spies * length(labeled_inds)))

  # build pseudo-labels \tilde{c}
  y1 <- y
  levels(y1) <- c(levels(y), paste0("U"))
  y1[c(unlabeled_inds, S)] <- "U"
  data[, toString(formula[[2]])] <- y1

  # train classifier
  mod <- BayesClassifier(formula, data, naive = naive, var_eps = var_eps)
  pred_S <- predict(mod, newdata = data[S,], type = "probs") # PREDICT SPIES
  t <- min(pred_S[cbind(1:length(S), y[S])]) # minimum probability of y in the correct class

  # train GMM on unlabeled data from new class
  pred_unlabeled <- predict(mod, newdata = data[unlabeled_inds,], type = "probs") # PREDICT UNLABELED
  newclass_inds <- unlabeled_inds[which(apply(pred_unlabeled[,setdiff(colnames(pred_unlabeled), "U")], 1, max) < t)] # PREDICTED CLASS HAS SMALL PROBABILITY

  # bic <- c()
  # models <- list()
  # for(g in 1:max_unknown){
  #   print(length(newclass_inds))
  #   try({
  gmm <- Mclust(X[newclass_inds,],
                  #G = g,
                  modelNames = ifelse(naive, "VVI", "VVV")) # train GMM on these new data
      print(paste(gmm$G, "new classes found"))

      # append new elements to classifier
      levels(y) <- c(levels(y), paste0("U",1:gmm$G))
      y[newclass_inds] <- paste0("U", predict(gmm, newdata = X[newclass_inds, ])$classification) # set new class elements
      data[, toString(formula[[2]])] <- y

      # re-train with all classes
      labeled_inds <- which(!is.na(y))
      mod <- BayesClassifier(formula, data[labeled_inds,], naive = naive, var_eps = var_eps)
    # })
#
#     bic <- c(bic, BIC(mod))
#     models <- append(models, list(mod))
#   }
#   print(bic)

  # return(models[[which.max(bic)]])
  return(mod)
}
