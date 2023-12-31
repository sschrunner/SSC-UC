% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSC.R
\name{EM}
\alias{EM}
\title{Train a semi-supervised Bayes classifier using EM algorithm}
\usage{
EM(
  formula,
  data,
  naive = FALSE,
  prior = "proportional",
  var_eps = 0.01,
  fixed_labels = NULL,
  maxiter = 20,
  verbose = TRUE
)
}
\arguments{
\item{formula}{a formula object explaining the functional relation between input variables and target}

\item{data}{a data.frame containing the dataset}

\item{naive}{whether a Bayes or Naive Bayes classifier should be used}

\item{prior}{type of class prior; either "uniform" (all classes are equally weighted) or "proportional" (all classes are weighted proportionally to their numbers of samples in the training data)}

\item{var_eps}{scalar to add to the main diagonal of the covariance matrix to assure numerical stability. If 0, no scalar will be added.}

\item{fixed_labels}{indices of rows in data, which have fixed labels that must not be changed during training}

\item{maxiter}{maximum number of iterations in EM algorithm}

\item{verbose}{whether outputs should be shown or suppressed}
}
\value{
a ´BayesClassifier´ object, see \link{BayesClassifier}
}
\description{
Builds and trains a semi-supervised classifier using the Expectation-Maximization (EM) algorithm.
}
\details{
Given a semi-supervised setup with labeled and unlabeled training samples, the function returns a BayesClassifier comprising all classes covered by the labeled training data.
All classes are assumed to be known and represented in the labeled training dataset.
Labeled training data are fixed (specified as ´fixed_labels´) and do not change during the algorithm. Unlabeled data are given with 'NA' as class label.
Training is performed using the Expectation-Maximization (EM) algorithm, which comprises two alternating steps: (a) estimating the model parameters of a BayesClassifier model from the current labeling, and (b) updating the labels by predicting from the current BayesClassifier model.
The algorithm terminates, if either a maximum number of iterations ´maxiter´ is reached, or if the same labels are predicted in two consecutive iterations.
}
\seealso{
\link{SSC} for semi-supervised classification with known and unknown classes, \link{BayesClassifier} for supervised classification
}
