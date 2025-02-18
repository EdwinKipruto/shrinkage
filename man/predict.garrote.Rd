% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonnegative_garrote.R
\name{predict.garrote}
\alias{predict.garrote}
\title{make predictions from a "garrote" object.}
\usage{
\method{predict}{garrote}(object, newx = NULL, ...)
}
\arguments{
\item{object}{Fitted \code{"garrote"} model object}

\item{newx}{Matrix of new values for \code{x} at which predictions are to be
made. Must be a matrix.}

\item{...}{not used at the moment}
}
\description{
Similar to other predict methods, this functions predicts fitted values,
from a fitted \code{"garrote"} object.
}
