#' Train an SSC-UC model
#' @description Trains a semi-supervised classifier with awareness of unknown classes, as described in \insertCite{Schrunner2020}{SSCUC}.
#' @details Given a semi-supervised setup with labeled and unlabeled training data, as well as known and unknown classes (classes represented in the labeled and unlabeled training data or in unlabeled training data only, respectively), a Bayes classifier shall be trained.
#' The algorithm is described in \insertCite{Schrunner2020}{SSCUC}.
#' @references
#'  \insertAllCited{}
#' @param formula a formula object explaining the functional relation between input variables and target
#' @param data a data.frame containing the dataset
#' @param perc_spies percentage of unlabeled points sampled as spies
#' @param naive whether a Bayes or Naive Bayes classifier should be used
#' @param prior type of class prior; either "uniform" (all classes are equally weighted) or "proportional" (all classes are weighted proportionally to their numbers of samples in the training data)
#' @param var_eps scalar to add to the main diagonal of the covariance matrix to assure numerical stability. If 0, no scalar will be added.
#' @param max_unknown maximum number of unknown classes (when number of unknown classes is determined automatically)
#' @param fixed_unknown fixed number of unknown classes (when number of unknown classes is specified manually)
#' @param runs number of bootstrap runs for selecting "likely unknowns"
#' @returns a ´BayesClassifier´ object, see \link{BayesClassifier}
#' @seealso \link{EM} for semi-supervised classification with known classes only, \link{BayesClassifier} for supervised classification
#' @importFrom Rdpack reprompt
#' @import mclust
#' @import stats
#' @export
SSC <- function(formula, data,
                perc_spies = 0.05, naive = FALSE, prior = "proportional",
                var_eps = 1e-2,
                max_unknown = 6, fixed_unknown = NULL,
                runs = 10){

  # TODO parameter consistency checks
  if(!prior %in% c("proportional", "uniform")){stop("Unknown prior option")}


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
    const_cols <- apply(X[newclass_inds,, drop = F], 2, function(x){return(length(unique(x)))}) == 1
    gmmModelName <- ifelse(naive, "VVI", "VVV")
    if(is.null(fixed_unknown)){
      g_opts <- 1:max_unknown
    }
    else{
      g_opts <- fixed_unknown
    }

    gmm <- Mclust(X[newclass_inds, !const_cols, drop = F],
                  G = g_opts,
                  modelNames = gmmModelName,
                  warn = TRUE) # train GMM on these new data

    if(is.null(gmm)){
      stop("Error: GMM did not run properly")
    }
    else{
      g <- gmm$G
      stop <- FALSE
      print(paste("Starting EM with", g, "classes"))

      # append new elements to classifier
      levels(y) <- c(known_classes, paste0("U",1:gmm$G))
      y[newclass_inds] <- paste0("U",
                                 predict(gmm, newdata = X[newclass_inds, !const_cols])$classification) # set new class elements
      data[, toString(formula[[2]])] <- y

      # re-train with all classes
      labeled_inds <- which(!is.na(y))
      mod <- EM(formula, data,
                naive = naive,
                prior = prior,
                var_eps = var_eps,
                fixed_labels = labeled_inds)
      print(paste("BIC:",BIC(mod)))

      # iterate to find better number of unknown classes, unless fixed by the user
      if(is.null(fixed_unknown)){
        while(!stop && g > 1){
          print(paste("Trying EM with", g-1, "classes"))
          gmm1 <- Mclust(X[newclass_inds, !const_cols, drop = F],
                         G = g-1,
                         modelNames = gmmModelName) # train alternative GMM
          if(is.null(gmm1)){
            warning("Error: GMM with fewer classes did not run properly - stopping")
            stop = TRUE
          }
          else{
            # append new elements to classifier
            y[newclass_inds] <- paste0("U",
                                       predict(gmm1, newdata = X[newclass_inds, !const_cols])$classification) # set new class elements
            y <- droplevels(y)
            data[, toString(formula[[2]])] <- y

            # re-train with all classes
            labeled_inds <- which(!is.na(y))
            mod1 <- EM(formula, data,
                       naive = naive,
                       prior = prior,
                       var_eps = var_eps,
                       fixed_labels = labeled_inds)
            print(paste("BIC:",BIC(mod1)))

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
  }

  return(mod)
}

#' Train a semi-supervised Bayes classifier using EM algorithm
#' @description Builds and trains a semi-supervised classifier using the Expectation-Maximization (EM) algorithm.
#' @details Given a semi-supervised setup with labeled and unlabeled training samples, the function returns a BayesClassifier comprising all classes covered by the labeled training data.
#' All classes are assumed to be known and represented in the labeled training dataset.
#' Labeled training data are fixed (specified as ´fixed_labels´) and do not change during the algorithm. Unlabeled data are given with 'NA' as class label.
#' Training is performed using the Expectation-Maximization (EM) algorithm, which comprises two alternating steps: (a) estimating the model parameters of a BayesClassifier model from the current labeling, and (b) updating the labels by predicting from the current BayesClassifier model.
#' The algorithm terminates, if either a maximum number of iterations ´maxiter´ is reached, or if the same labels are predicted in two consecutive iterations.
#' @inheritParams SSC
#' @param fixed_labels indices of rows in data, which have fixed labels that must not be changed during training
#' @param maxiter maximum number of iterations in EM algorithm
#' @param verbose whether outputs should be shown or suppressed
#' @returns a ´BayesClassifier´ object, see \link{BayesClassifier}
#' @seealso \link{SSC} for semi-supervised classification with known and unknown classes, \link{BayesClassifier} for supervised classification
#' @export
EM <- function(formula, data,
               naive = FALSE,
               prior = "proportional",
               var_eps = 1e-2,
               fixed_labels = NULL,
               maxiter = 20,
               verbose = TRUE){

  editable_labels <- setdiff(1:nrow(data), fixed_labels)
  if(length(editable_labels) == 0){
    stop("Not all labels can be fixed")
  }

  stop <- FALSE
  iter <- 1

  while(!stop && iter < maxiter){
    y <- convertFormula(formula, data)$y

    # E-step (with non-NA data only)
    model <- BayesClassifier(formula, data,
                             naive = naive,
                             prior = prior,
                             var_eps = var_eps)

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
