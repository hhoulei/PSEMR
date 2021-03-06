% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSE_IVW_order.R
\name{Ordered_IVW}
\alias{Ordered_IVW}
\title{Causal Mediation Analysis with Multiple Causally Ordered Mediators}
\usage{
Ordered_IVW(betaYG, sebetaYG, betaXG, sebetaXG, betaMG, sebetaMG)
}
\arguments{
\item{betaYG}{-- a vector, genetic summary statistics for outcome (coefficient)}

\item{sebetaYG}{-- a vector, genetic summary statistics for outcome (standard error)}

\item{betaXG}{-- a vector, genetic summary statistics for exposure (coefficient)}

\item{sebetaXG}{-- a vector,genetic summary statistics for exposure (standard error)}

\item{betaMG}{-- matrix or dataframe with n mediators, genetic summary statistics for mediators (coefficient)}

\item{sebetaMG}{-- matrix or dataframe with n mediators, genetic summary statistics formediators (standard error)}
}
\value{
list
}
\description{
Causal Mediation Analysis with Multiple Causally Ordered Mediators
}
