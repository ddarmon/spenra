% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/er-functions.R
\name{timeseries.entropyrate.lwo}
\alias{timeseries.entropyrate.lwo}
\title{Estimate the Entropy Rate of a Scalar Time Series Using a Plug-In Leave-Window-Out Estimator}
\usage{
timeseries.entropyrate.lwo(x, pred.dens, half.window.length = 1)
}
\arguments{
\item{x}{The scalar time series.}

\item{pred.dens}{The npcdens-based predictive density used to estimate the specific entropy rate.}

\item{half.window.length}{The half-length of the window about the present removed in estimating the specific entropy rate.}
}
\description{
timeseries.entropyrate.lwo estimates the specific entropy rate using a plug-in, leave-window-out estimator.
}

