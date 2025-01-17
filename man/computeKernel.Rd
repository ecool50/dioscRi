% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMD_Vae.R
\name{computeKernel}
\alias{computeKernel}
\title{Compute RBF Kernel for MMD Calculation}
\usage{
computeKernel(x, y)
}
\arguments{
\item{x}{Tensor; first input tensor.}

\item{y}{Tensor; second input tensor.}
}
\value{
Tensor containing the RBF kernel matrix.
}
\description{
Computes the RBF kernel between two inputs, typically for use in MMD calculations.
}
