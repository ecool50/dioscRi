% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling_and_visualization.R
\name{plotAUC}
\alias{plotAUC}
\title{Plot AUC for Model Evaluation}
\usage{
plotAUC(fit, xTest, yTest, title = "")
}
\arguments{
\item{fit}{Fitted model object.}

\item{xTest}{Matrix of test data.}

\item{yTest}{Vector of true labels.}

\item{title}{Character; plot title.}
}
\value{
A list with a ggplot object for the ROC curve and a data frame of predictions.
}
\description{
Plots the ROC curve and calculates AUC for model predictions. Also returns
a data frame with prediction results.
}
