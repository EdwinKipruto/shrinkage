% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lasso.R
\name{glmnetCoef}
\alias{glmnetCoef}
\title{Fits ridge or lasso regression model and the resulting estimates are used by garrote and
adaptive lasso as initial estimates}
\usage{
glmnetCoef(
  x,
  y,
  nfolds = 10,
  foldid = NULL,
  type.measure = c("mse", "mae"),
  alpha = 1,
  lambda = NULL,
  nlambda = 100,
  standardize = FALSE,
  lambda.min.ratio = 1e-04
)
}
\arguments{
\item{x}{A standardized matrix of predictors of dimension nobs x nvars, where
each row is an observation vector.If not standardized, set standardize to TRUE.}

\item{y}{A numeric centered quantitative vector of response. If not centered,
set `standardize` to TRUE; this will standardize x but center y. See `standardize.response`
for more details.}

\item{nfolds}{Number of folds for cross-validation. Default is 10. Smallest value
allowable is nfolds = 3.}

\item{foldid}{an optional vector of values between 1 and nfolds identifying
what fold each observation is in. If supplied, nfolds can be missing.}

\item{type.measure}{Loss function to use for cross-validation. Currently,
two options are available. The default option is \code{type.measure = "mse"},
which corresponds to squared error. Another option is \code{type.measure = "mae"}
(mean absolute error).}

\item{alpha}{The alpha parameter controls the mixing of ell1 and ell2 penalties.
It is a value between 0 and 1, where 0 represents ridge regression and 1 represents lasso
regression.}

\item{lambda}{Optional user-supplied lambda sequence. Default is NULL, and
the program chooses its own sequence assuming x and y are standardized.}

\item{nlambda}{The number of lambda values.Default is 100}

\item{standardize}{Specifies whether to standardize the predictors matrix
\code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
standardized. Additionally, the response variable \code{`y`} will be centered
by subtracting its mean. This standardization step can be useful to ensure
that the predictors are on a comparable scale. By default,
\code{standardize = FALSE}, indicating that no standardization will be performed.
This assumes that users have already standardized x and centered y so that the intercept
is zero.}

\item{lambda.min.ratio}{When lambda values are automatically generated, the sequence
is determined by lambda.max and lambda.min ratio. The program generates nlambda values
on the log scale from lambda.max down to lambda.min. lambda.max is not user-specified but is
computed from the input standardized x and y. It is the smallest value for lambda such that
all the coefficients are zero. The default is lambda.min.ratio = 0.0001.}
}
\value{
A list containing the following items:
\item{beta:}{The estimated regression coefficients of the model without an intercept.}
\item{lambda.initial:}{The tuning parameters used to fit the model.}
}
\description{
Fits ridge or lasso models with the optimal tuning parameter determined
through cross-validation. If the number of nonzero coefficients is less than
two, the optimal lambda is replaced with a lambda value that results in at
least two nonzero coefficients. This adjustment is necessary because the
nonnegative garrote and adaptive lasso programs use the glmnet function,
which requires at least two predictor variables in the model.
}
