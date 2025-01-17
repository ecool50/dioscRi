% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling_and_visualization.R
\name{selectLambda}
\alias{selectLambda}
\title{Select Optimal Lambda for Group Lasso Model}
\usage{
selectLambda(
  xTrain,
  yTrain,
  groups,
  penalty = "cMCP",
  alphaRange = seq(0.1, 1, length = 10),
  lambdaSearch = seq(0.01, 1, length = 10),
  seed = 1994
)
}
\arguments{
\item{xTrain}{Matrix of training features.}

\item{yTrain}{Vector of training labels.}

\item{groups}{List of groups for overlapping group lasso.}

\item{penalty}{Character; penalty type, default is "cMCP".}

\item{lambdaSearch}{Numeric vector of lambda values to search.}

\item{seed}{Integer; random seed for reproducibility.}

\item{weights}{Vector of weights for the lasso penalty.}

\item{alpha}{Numeric; fixed alpha value for the group lasso.}
}
\value{
A list with the best fitted model, best BIC, best lambda, and BIC values for each lambda.
}
\description{
Searches for the best lambda value for the group lasso model by minimizing BIC.
}
