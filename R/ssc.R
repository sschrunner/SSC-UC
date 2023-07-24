#' Train an SSC-UC model
#' @description Builds and trains a semi-supervised classifier with awareness of unknown classes, as described in Schrunner et al. (2020).
#' @param formula a formula object explaining the functional relation between input variables and target
#' @param data a data.frame containing the dataset
#' @param perc_spies percentage of unlabeled points sampled as spies
#' @param naive whether a Bayes or Naive Bayes classifier should be used
#' @param var_eps scalar to add to the main diagonal of the covariance matrix to assure numerical stability. If 0, no scalar will be added.
#' @param max_unknown maximum number of unknown classes
#' @param runs number of bootstrap runs for selecting "likely unknowns"
#' @return a BayesClassifier object
#' @import mclust
#' @import stats
#' @export
SSC <- function(formula, data,
                perc_spies = 0.05, naive = FALSE, var_eps = 1e-2, max_unknown = 6,
                runs = 10){
  l <- convertFormula(formula, data)
  X <- l$X
  y <- l$y
  known_classes <- levels(y)

  # extract information on labeled and unlabeled samples
  labeled_inds <- which(!is.na(y))
  num_labeled <- length(labeled_inds)
  unlabeled_inds <- which(is.na(y))
  num_unlabeled <- length(unlabeled_inds)

  spy_selection <- function(){
    # sample spies
    S <- sample(labeled_inds, size = round(perc_spies * length(labeled_inds)))

    # build pseudo-labels \tilde{c}
    y1 <- y
    levels(y1) <- c(known_classes, paste0("U"))
    y1[c(unlabeled_inds, S)] <- "U"
    data[, toString(formula[[2]])] <- y1

    # train classifier
    mod <- BayesClassifier(formula, data, naive = naive, var_eps = var_eps)
    pred_S <- predict(mod, newdata = data[S,], type = "probs") # PREDICT SPIES
    t <- quantile(pred_S[cbind(1:length(S), y[S])], probs = 0.01) # minimum or low quantile of probabilities of y in the correct class

    # train GMM on unlabeled data from new class
    pred_unlabeled <- predict(mod, newdata = data[unlabeled_inds,], type = "probs") # PREDICT UNLABELED
    newclass_inds <- unlabeled_inds[which(apply(pred_unlabeled[,setdiff(colnames(pred_unlabeled), "U")], 1, max) < t)] # PREDICTED CLASS HAS SMALL PROBABILITY
    return(newclass_inds)
  }

  # multiple bootstrap-runs over spy_selection
  bootstrap <- replicate(runs, spy_selection())
  # union over all likely unknowns
  newclass_inds <- unlist(bootstrap)
  # apply 50% threshold across bootstrap samples
  newclass_inds <- unique(newclass_inds)[table(newclass_inds) > (runs / 2)]

  # only if new classes were detected
  if(length(newclass_inds) > 0){
    # bic <- c()
    # models <- list()
    # for(g in 1:max_unknown){
    #   print(length(newclass_inds))
    #   try({
    uniform_cols <- apply(X[newclass_inds,, drop = F], 2, function(x){return(length(unique(x)))}) == 1
    gmmModelName <- ifelse(naive, "VVI", "VVV")

    gmm <- Mclust(X[newclass_inds, !uniform_cols, drop = F],
                  modelNames = gmmModelName) # train GMM on these new data

    if(is.null(gmm)){
      stop("Error: GMM did not run properly")
    }
    else{
      g <- gmm$G
      stop <- FALSE
      print(paste("Starting EM with", g, "classes"))

      # append new elements to classifier
      levels(y) <- c(known_classes, paste0("U",1:gmm$G))
      y[newclass_inds] <- paste0("U", predict(gmm, newdata = X[newclass_inds, ])$classification) # set new class elements
      data[, toString(formula[[2]])] <- y

      # re-train with all classes
      labeled_inds <- which(!is.na(y))
      #mod <- BayesClassifier(formula, data[labeled_inds,], naive = naive, var_eps = var_eps)
      mod <- EM(formula, data,
                naive = naive, var_eps = var_eps,
                fixed_labels = labeled_inds)
      print(paste("BIC:",BIC(mod)))

      while(!stop && g > 1){
        print(paste("Trying EM with", g-1, "classes"))
        gmm1 <- Mclust(X[newclass_inds, !uniform_cols, drop = F],
                       G = g-1,
                       modelNames = gmmModelName) # train alternative GMM
        if(is.null(gmm1)){
          warning("Error: GMM with fewer classes did not run properly - stopping")
          stop = TRUE
        }
        else{
          # append new elements to classifier
          y[newclass_inds] <- paste0("U", predict(gmm1, newdata = X[newclass_inds, ])$classification) # set new class elements
          y <- droplevels(y)
          data[, toString(formula[[2]])] <- y

          # re-train with all classes
          labeled_inds <- which(!is.na(y))
          #mod1 <- BayesClassifier(formula, data[labeled_inds,], naive = naive, var_eps = var_eps)
          mod1 <- EM(formula, data,
                     naive = naive, var_eps = var_eps,
                     fixed_labels = labeled_inds)
          print(paste("BIC:",BIC(mod1)))
          #return(list(mod =mod, mod1 = mod1))

          if(BIC(mod1) >= BIC(mod)){
            print("BIC increased when updating - stopping")
            stop <- TRUE
          }
          else{
            mod <- mod1
            g <- g-1
          }
        }
      }
    }
  }

  return(mod)
}

#' Train a semi-supervised Bayes classifier using EM
#' @description Builds and trains a semi-supervised classifier using the Expectation-Maximization (EM) algorithm.
#' @inheritParams SSC
#' @param fixed_labels indices of rows in data, which have fixed labels that must not be changed during training
#' @param verbose whether outputs should be shown or suppressed
#' @return a list containing a BayesClassifier object and a vector of labels
#' @export
EM <- function(formula, data,
               naive = FALSE, var_eps = 1e-2,
               fixed_labels = NULL,
               verbose = TRUE){

  editable_labels <- setdiff(1:nrow(data), fixed_labels)
  if(length(editable_labels) == 0){
    stop("Not all labels can be fixed")
  }

  stop <- FALSE
  iter <- 1

  while(!stop && iter < 10){
    y <- convertFormula(formula, data)$y

    # E-step
    model <- BayesClassifier(formula, data,
                             naive = naive, var_eps = var_eps)
    # M-step
    y[editable_labels] <- predict(model,
                                  newdata = data[editable_labels,],
                                  type = "class")

    # if new class labels are equal to old ones
    if(iter > 1 && all(y == data[, toString(formula[[2]])])){
      stop <- TRUE
    }
    else{
      data[, toString(formula[[2]])] <- y
      iter <- iter + 1
    }
  }

  if(verbose){print(paste("EM converged after", iter, "iterations"))}
  return(model)
}
