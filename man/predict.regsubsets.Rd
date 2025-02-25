% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_selection.R
\name{predict.regsubsets}
\alias{predict.regsubsets}
\title{Predict Method for regsubsets}
\usage{
\method{predict}{regsubsets}(object, newx = NULL, id, ...)
}
\arguments{
\item{object}{Fitted "regsubsets" model.}

\item{newx}{Matrix of new values for 'x' at which predictions are to be made. Must be a matrix.}

\item{id}{The model size at which predictions are required. Note that different model sizes are returned by regsubsets.}

\item{...}{not used at the moment}
}
\description{
Make predictions from a "regsubsets" object. Similar to other predict methods.
}
