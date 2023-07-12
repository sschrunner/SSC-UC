#' Train an SSC-UC model
#' @description Builds and trains a semi-supervised classifier with awareness of unknown classes, as described in Schrunner et al. (2020).
#' @param formula a formula object explaining the functional relation between input variables and target
#' @param data a data.frame containing the dataset
#' @param perc_spies percentage of unlabeled points sampled as spies
#' @param naive whether a Bayes or Naive Bayes classifier should be used
#' @param var_eps scalar to add to the main diagonal of the covariance matrix to assure numerical stability
#' @return a BayesClassifier object
#' @import mclust
#' @import stats
#' @export
SSC <- function(formula, data,
                perc_spies = 0.05, naive = FALSE, var_eps = 1e-2){
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y

  # extract information on labeled and unlabeled samples
  labeled_inds <- which(!is.na(y))
  num_labeled <- length(labeled_inds)
  unlabeled_inds <- which(is.na(y))
  num_unlabeled <- length(unlabeled_inds)

  # class labels
  N <- length(levels(y))
  prior <- table(y) / length(labeled_inds)

  # implement training function
  S <- sample(labeled_inds, size = round(perc_spies * length(labeled_inds)))

  # build pseudo-labels \tilde{c}
  y1 <- y
  levels(y1) <- c(levels(y), paste0("U"))
  y1[c(unlabeled_inds, S)] <- "U"
  data[, toString(formula[[2]])] <- y1

  # train classifier
  mod <- BayesClassifier(formula, data, naive = naive)
  pred_S <- predict(mod, newdata = data[S,], type = "probs") # PREDICT SPIES
  t <- min(pred_S[cbind(1:length(S), y[S])]) # minimum probability of y in the correct class

  # train GMM on unlabeled data from new class
  pred_unlabeled <- predict(mod, newdata = data[unlabeled_inds,]) # PREDICT UNLABELED
  newclass_inds <- unlabeled_inds[which(apply(pred_unlabeled, 1, max) < t)] # PREDICTED CLASS HAS SMALL PROBABILITY
  gmm <- Mclust(X[newclass_inds,],
                G = 1:10,
                modelNames = ifelse(naive, "VVI", "VVV")) # train GMM on these new data
  # TODO: add BIC to select best GMM model

  # append new elements to classifier
  prior <- c(U = sum(pred_S < t) / length(pred_S), prior) # prior for unknown class
  levels(y) <- c(levels(y), paste0("U",1:gmm$G))
  y[newclass_inds] <- paste0("U", predict(gmm, newdata = X[newclass_inds, ])$classification) # set new class elements
  data[, toString(formula[[2]])] <- y
  print(paste(gmm$G, "new classes found"))

  # re-train with all classes
  labeled_inds <- which(!is.na(y))
  mod <- BayesClassifier(formula, data[labeled_inds,], naive = naive)

  return(mod)
}
