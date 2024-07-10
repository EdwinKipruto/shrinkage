#' @title Fit Gaussian linear model with lasso regularization
#'
#' @description Fits the Lasso model using the optimal tuning parameter
#' determined through cross-validation. It also supports the selection
#' of tuning parameters based on the Akaike Information Criterion (AIC) or
#' Bayesian Information Criterion (BIC). The value of lambda value that yield the
#' smallest AIC or BIC is chosen as the best tuning parameter.
#' @param x A standardized matrix of predictors of dimension nobs x nvars, where
#' each row is an observation vector.If not standardized, set standardize to TRUE.
#' @param y A numeric centered quantitative vector of response. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. See `standardize.response`
#' for more details.
#' @param standardize Specifies whether to standardize the predictors matrix
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#' standardized. Additionally, the response variable \code{`y`} will be centered
#' by subtracting its mean. This standardization step can be useful to ensure
#' that the predictors are on a comparable scale. By default,
#' \code{standardize = FALSE}, indicating that no standardization will be performed.
#' This assumes that users have already standardized x and centered y so that the intercept
#' is zero.
#' @param standardize.response Specifies whether the response variable y should be
#' standardized to have unit variance. This option divides \code{`y`} by its standard deviation.
#' @param foldid an optional vector of values between 1 and nfolds identifying
#' what fold each observation is in. If supplied, nfolds can be missing.
#' @param nfolds Number of folds for cross-validation. Default is 10. Smallest value
#' allowable is nfolds = 3.
#' @param lambda Optional user-supplied lambda sequence. Default is NULL, and
#' the program chooses its own sequence assuming x and y are standardized.
#' @param nlambda The number of lambda values.Default is 100
#' @param type.measure Loss function to use for cross-validation. Currently,
#' two options are available. The default option is \code{type.measure = "mse"},
#' which corresponds to squared error. Another option is \code{type.measure = "mae"}
#' (mean absolute error).
#' @param criterion The criterion used to select tuning parameters for the
#' lasso regression. Available options include: \code{criterion = "cv"}
#' which perform cross-validation to select the optimal tuning parameter, while
#' \code{criterion = "aic"} and \code{criterion = "bic"} select the tuning parameter
#'  based on the AIC and BIC, respectively.
#' @param sigma2 The residual variance obtained by fitting the full OLS model without
#' any regularization or feature selection. It is used for the calculation of AIC and BIC.
#' To compute these metrics, the residual variance needs to be provided via the
#' `sigma2` parameter.
#' @param parallel The program use \code{"[cv.glmnet()]"} function which supports
#' parallel computing. To make it work, users must register parallel beforehand.
#' see `glmnet` package for details
#' @param lambda.min.ratio When lambda values are automatically generated, the sequence
#' is determined by lambda.max and lambda.min ratio. The program generates nlambda values
#' on the log scale from lambda.max down to lambda.min. lambda.max is not user-specified but is
#' computed from the input standardized x and y. It is the smallest value for lambda such that
#' all the coefficients are zero. The default is lambda.min.ratio = 0.0001.
#' @param betatypes Not used but added for consistency with the oracle_model().
#' @return Returns the following items:
#' \item{beta:}{Shrunken regression estimates without an intercept.}
#' \item{lambda:}{The tuning parameter used for the estimation of parameters.}
#' \item{gamma:}{Not applicable for lasso but returned for consistency with other methods like adaptive lasso.}
#' \item{x:}{A standardized matrix of predictors used in fitting the linear model.}
#' \item{y:}{A centered vector of the response variable used in fitting the linear model.}
#' @import glmnet
#' @export
lassofit <- function(x, y,
                     nfolds = 10,
                     foldid = NULL,
                     lambda = NULL,
                     nlambda = 100,
                     type.measure = c("mse", "mae"),
                     parallel = FALSE,
                     criterion = c("cv", "aic", "bic"),
                     standardize = FALSE,
                     standardize.response = FALSE,
                     sigma2 = NULL,
                     lambda.min.ratio = 0.0001,
                     betatypes = NULL){

  # match arguments and set up defaults
  type.measure <- match.arg(type.measure)
  criterion <- match.arg(criterion)

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # x must have column names
  if (is.null(colnames(x)))
    stop("x must have column names.")

  # assert that residual variance must be provided for AIC and BIC criteria
  if (criterion != "cv" && is.null(sigma2)) {
    stop(sprintf("!sigma2 must be provided for %s calculation.", toupper(criterion)))
  }

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)

    # center x and scale x by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }

  # Standardize the response variable to have sd of 1
  if(standardize.response){
    y <- y / sd_x(y)
  }

  # Generate sequence of lambdas when not provided.
  if (is.null(lambda)) {
    lambda <- lambda_seq_lasso(x = x,
                               y = y,
                               nlambda = nlambda,
                               epsilon = lambda.min.ratio
                              )
  }

  # override the default nlambda when lambdas are supplied by the user
  nlambda <- length(lambda)

  # regression estimates without intercept for selected criterion. We assume
  # x is standardized and y centered
  fit <- switch (criterion,
    "cv" = lasso_cv(x = x, y = y,
                    nfolds = nfolds,
                    foldid = foldid,
                    lambda = lambda,
                    nlambda = nlambda,
                    parallel = parallel,
                    standardize = FALSE,
                    type.measure =  type.measure
                   ),

    "aic" = lasso_aic_bic(x = x, y = y,
                          lambda = lambda,
                          nlambda = nlambda,
                          standardize = F,
                          criterion = "aic",
                          sigma2 = sigma2
                         ),

    "bic" = lasso_aic_bic(x = x, y = y,
                          lambda = lambda,
                          nlambda = nlambda,
                          standardize = FALSE,
                          criterion = "bic",
                          sigma2 = sigma2)
  )

  # return regression estimates, x and y used to fit the model
  fit_parms <- list(beta = fit$beta,
                    lambda =fit$lambda,
                    gamma = NA,
                    x = x, y = y)

  # set class to be used by generic functions
  class(fit_parms) = "lassofit"

  return(fit_parms)
}

#-----Lasso with cross-validation-----------------------------------------------

#' Helper function that fits lasso with cross-validation tuning parameters
#' @inheritParams lassofit
lasso_cv <- function(x, y,
                     nfolds,
                     foldid,
                     lambda,
                     nlambda,
                     type.measure,
                     parallel,
                     standardize){

  # sample size
  n <- dim(x)[1]

  # set foldid for splitting data in cross-validation
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # Estimate optimal tuning parameter and estimate corresponding parameters.
  optimal_lambda <- cv.glmnet(x = x, y = y,
                              foldid = foldid,
                              family = "gaussian",
                              lambda = lambda,
                              alpha= 1,
                              nlambda = nlambda,
                              type.measure = type.measure,
                              parallel = parallel,
                              standardize = standardize)$lambda.min

  fit <- glmnet(x = x, y = y,
                lambda = optimal_lambda,
                family = "gaussian",
                standardize = standardize,
                alpha = 1)

  # regression estimates without the intercept
  out <- list(beta = setNames(as.numeric(fit$beta), colnames(x)),
              lambda =optimal_lambda)
  return(out)
  }

#-----Lasso with AIC and BIC tuning parameters----------------------------------

#' Helper function that fits lasso with aic and bic tuning parameters
#' @inheritParams lassofit
lasso_aic_bic <- function(x, y,
                          lambda,
                          nlambda,
                          standardize,
                          criterion,
                          sigma2){

  # sample size and number of variables
  np <- dim(x)
  n <- as.integer(np[1L])
  p <- as.integer(np[2L])

  # Note that glmnet supports vector of lambdas
  fit <- glmnet(x = x, y = y,
                lambda = lambda,
                standardize = standardize,
                alpha = 1,
                family = "gaussian"
               )

  # Matrix of regression coefficients without intercept
  betas <- fit$beta

  # linear predictor: crossprod in Matrix package supports sparse matrices
  xb <- as.matrix(Matrix::crossprod(t(x), betas))

  # residual for each value of lambda: each column is for a specific lambda
  res <- matrix(y, nrow = nrow(xb), ncol = ncol(xb), byrow = FALSE) - xb

  # a vector of residual sum of squares for each value of lambda
  rss <- colSums(res ^ 2, na.rm = TRUE)

  # degrees of freedom without intercept: use rank(xa) instead of number of
  # nonzero coefficients, applicable for general x matrix. It can reduce nonzero
  # coefficients when x is full rank or in low correlation
  # see degrees of freedom in lasso problems by Tibshirani and Jonathan (2012)
  # DOI: 10.1214/12-AOS1003
  #df <- fit$df

  # vector of df for each value of lambda
  df <- calculate_df(x = x, betas = as.matrix(betas))

  # compute AIC and BIC for each value of lambda
  pn <- (rss / (n * sigma2))
  aic <- pn + (2 / n) * df
  bic <- pn + (log(n) / n) * df

  # Find lambda that yield smallest AIC and BIC
  lambda.aic <- lambda[which.min(aic)]
  lambda.bic <- lambda[which.min(bic)]

  # beta <- switch(criterion,
  #              "aic" = as.numeric(predict(fit, type = "coefficients",
  #                                         s = lambda.aic,
  #                                         exact = TRUE,
  #                                         alpha = 1,
  #                                         standardize = standardize,
  #                                         x = x,
  #                                         y = y))[-1],
  #              "bic" = as.numeric(predict(fit, type = "coefficients",
  #                                         s = lambda.bic,
  #                                         exact = TRUE,
  #                                         alpha = 1,
  #                                         standardize = standardize,
  #                                         x = x,
  #                                         y = y))[-1]
  #   )

  # Similar to the above code but no warnings
  beta <- switch(criterion,
                 "aic" = as.numeric(glmnet(x = x, y = y,
                                           lambda = lambda.aic,
                                           standardize = standardize,
                                           alpha = 1,
                                           family = "gaussian")$beta),
                 "bic" = as.numeric(glmnet(x = x, y = y,
                                           lambda = lambda.bic,
                                           standardize = standardize,
                                           alpha = 1,
                                           family = "gaussian")$beta)
  )

  fits <- list(beta = setNames(beta, colnames(x)),
    lambda = ifelse(criterion == "aic",
                    lambda.aic,
                    lambda.bic)
  )

   return(fits)
}

# Function that calculates lambda sequence for the lasso see
# https://stackoverflow.com/questions/25257780/how-does-glmnet-compute-the-maximal-lambda-value/

#' calculates sequence of lambdas for the lasso
#' @param x a standardized matrix of predictors
#' @param y a centered vector of response
#' @param nlambda Number of lambdas
#' @param epsilon The ratio of smallest value of the generated lambda
#' sequence (say lambda.min) to lambda.max. See lambda.min.ratio
#' in \code{[cv.glmnet]} for more details
#' @export
lambda_seq_lasso <- function(x, y, nlambda = 100, epsilon =  0.0001){
  # calculate lambda-max that gives all estimates equal to 0
  #lambda_max <- max(abs(t(x)%*%y))/length(y)
  #lambda_max <- max(abs(colSums(as.matrix(x) * y, na.rm = T))) / length(y)
  x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  lambda_max <- max(abs(Matrix::crossprod(x, y)), na.rm = TRUE) / length(y)

  # calculate lambda_min that gives all estimates equal to nonzero
  lambda_min <- epsilon * lambda_max

  # choose the grid of lambdas to be equidistant on the log-scale
  lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

  # return lambdas
  return(lambda)
}

#' Fits ridge or lasso regression model and the resulting estimates are used by garrote and
#' adaptive lasso as initial estimates
#'
#' @description Fits ridge or lasso models with the optimal tuning parameter determined
#' through cross-validation. If the number of nonzero coefficients is less than
#' two, the optimal lambda is replaced with a lambda value that results in at
#' least two nonzero coefficients. This adjustment is necessary because the
#' nonnegative garrote and adaptive lasso programs use the glmnet function,
#' which requires at least two predictor variables in the model.
#' @inheritParams lassofit
#' @param alpha The alpha parameter controls the mixing of ell1 and ell2 penalties.
#' It is a value between 0 and 1, where 0 represents ridge regression and 1 represents lasso
#' regression.
#' @returns A list containing the following items:
#' \item{beta:}{The estimated regression coefficients of the model without an intercept.}
#' \item{lambda.initial:}{The tuning parameters used to fit the model.}
#' @export
glmnetCoef<-function(x, y,
                     nfolds = 10,
                     foldid = NULL,
                     type.measure = c("mse", "mae"),
                     alpha = 1,
                     lambda = NULL,
                     nlambda = 100,
                     standardize = FALSE,
                     lambda.min.ratio = 0.0001){

  # match arguments and set up defaults
  type.measure <- match.arg(type.measure)

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # sample size
  dimx <- dim(x)
  n <- as.integer(dimx[1L])
  nv <- as.integer(dimx[2L])

  # check whethers dimensions are identical
  ny <- length(y)
  if (n!=ny) {
    stop(sprintf("The length of y (%d) is not equal to the number of rows of x (%d)",ny, n))
  }
  # set foldid for splitting data in cross-validation
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)

    # center x by mean and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }

  # assert that alpha must be either 0 (ridge) or 1 (lasso)
  if (!alpha %in% c(0, 1))
    stop("The value of alpha must be either 0 or 1. Other values are not
         supported at the moment.")

  # generate sequence of lambdas if not provided
  if (is.null(lambda)) {
    if (alpha == 1){

      lambda <- lambda_seq_lasso(x = x, y = y,
                                 nlambda = nlambda,
                                 epsilon =  lambda.min.ratio
                                )
    } else {
      lambda <- lambda_seq_ridge(x = x, y = y,
                                 nlambda = nlambda,
                                 epsilon =  lambda.min.ratio
                                )
    }
  }

  # length of lambdas
  nlambda <- length(lambda)

  # lasso or ridge optimal solutions. if x is standardized and y centered the
  # estimated intercept will be approx 0
  lambda_min <- lambda

  if (length(lambda)>1){
  lambda_min <- glmnet::cv.glmnet(x = x, y = y,
                                  type.measure = type.measure,
                                  alpha = alpha, foldid =foldid,
                                  lambda = lambda, family = "gaussian",
                                  nlambda = nlambda,weights = NULL,
                                  offset = NULL,
                                  standardize = FALSE)$lambda.min
  }

  # Glmnet fit
  # fit <- glmnet(x = x, y = y,
  #               alpha = alpha,
  #               family = "gaussian",
  #               lambda = lambda,
  #               weights = NULL,
  #               offset = NULL,
  #               standardize = FALSE)

  # Optimal regression coefficients without intercept
  # beta <- as.numeric(predict(fit,
  #                            type = "coefficients",
  #                            x = x,
  #                            y = y,
  #                            s = lambda_min,
  #                            exact = T,
  #                            alpha = alpha,
  #                            family = "gaussian",
  #                            weights = NULL,
  #                            offset = NULL,
  #                           standardize = FALSE))[-1]

  # I prefer this option to have identical results with usual fit of lasso without using predict
  beta <- as.numeric(glmnet(x = x, y = y,
                alpha = alpha,
                family = "gaussian",
                lambda = lambda_min,
                weights = NULL,
                offset = NULL,
                standardize = FALSE)$beta)

 # Check if the number of nonzero coefficients is less than 2. If so, search for
 # a lambda value that produces exactly 2 nonzero coefficients. This is
 # particularly important for methods such as garrote and adaptive lasso since
 # it uses the resulting estimates as initial estimates. glmnet will then be used
 # but requires at least 2 variables to work

  if(sum(beta!=0)<2){
     fit <- glmnet::glmnet(x = x, y = y,
                           alpha = alpha,
                           family = "gaussian",
                           lambda = lambda, # we assume lambda is a vector
                           weights = NULL,
                           offset = NULL,
                           standardize = FALSE)
    # lambda that yields 2 or more nonzero coefficients
    lambda_min <- lambda[min(which(fit$df >= 2))]

    # beta <- as.numeric(predict(fit,
    #                            type = "coefficients",
    #                            s = lambda_min,
    #                            exact = TRUE,
    #                            x = x,
    #                            y = y,
    #                            weights = NULL,
    #                            offset = NULL,
    #                            alpha = alpha,
    #                            family = "gaussian",
    #                            standardize = FALSE))[-1]

    beta <- as.numeric(glmnet(x = x, y = y,
                              alpha = alpha,
                              family = "gaussian",
                              lambda = lambda_min,
                              weights = NULL,
                              offset = NULL,
                              standardize = FALSE)$beta)
  }

  out <- list(coefficients = setNames(beta, colnames(x)),
              lambda.initial = lambda_min)

  return(out)
}

#' Extract coefficients from a lassofit object
#'
#'@method coef lassofit
#'@param object Fitted \code{"lassofit"} model object
#'@param ... not used at the moment
#'@return regression estimates without intercept
#'@export
coef.lassofit <- function(object,...) {
  object$beta
}

#' make predictions from a "lassofit" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"lassofit"} object.
#' @param object Fitted \code{"lassofit"} model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict lassofit
#' @export
predict.lassofit <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  #newx = cbind(rep(1,dim(newx)[1L]),newx)
  newx %*% coef.lassofit(object)
}
