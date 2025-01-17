% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_processing.R
\name{trainCellTypeClassifier}
\alias{trainCellTypeClassifier}
\title{Train Cell Type Classifier using Caret}
\usage{
trainCellTypeClassifier(trainX, testX, model = "lda")
}
\arguments{
\item{trainX}{Data frame of training features with cell type labels.}

\item{testX}{Data frame of test features.}

\item{model}{Character; the method to use for training, e.g., "lda".}
}
\value{
Predicted cell types for the test set.
}
\description{
Trains a cell type classification model using LDA (or specified model) and
evaluates it on test data. Uses cross-validation for training.
}
