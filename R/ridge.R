#' @title Fit a linear model with ridge regularization
#'
#' @description Fits the ridge model using the optimal tuning parameter
#' determined through cross-validation. It also supports the selection
#' of tuning parameters based on the Akaike Information Criterion (AIC) or
#' Bayesian Information Criterion (BIC). The value of the lambda parameter
#' that yields the smallest AIC or BIC is chosen as the best tuning parameter.
#' @param x A standardized matrix of predictors of dimension nobs x nvars;
#' each row is an observation vector.If not standardized, set standardize to
#' TRUE.
#' @param y A numeric, centered quantitative vector of responses. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. See `standardize.response`
#' for more details.
#' @param standardize Specifies whether to standardize the predictors matrix
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. This standardization step can be useful to ensure
#'  that the predictors are on a comparable scale. By default,
#'  \code{standardize = FALSE}, indicating that no standardization will be performed.
#'  This assumes that users have already standardized x and centered y so that the intercept
#'  is zero.
#' @param standardize.response Specifies whether the response variable y should be
#' standardized to have unit variance. This option divides \code{`y`} by its standard
#' deviation.
#' @param foldid an optional vector of values between 1 and nfolds identifying
#' what fold each observation is in. If supplied, nfolds can be missing.
#' @param nfolds Number of folds for cross-validation. Default is 10. Smallest
#' value allowable is nfolds = 3.
#' @param lambda Optional user-supplied lambda sequence. Default is NULL, and
#' the program chooses its own sequence.
#' @param nlambda The number of lambda values. Default is 100
#' @param type.measure Loss function to use for cross-validation. Currently,
#' two options are available. The default option is \code{type.measure = "mse"},
#' which corresponds to squared error. Another option is \code{type.measure = "mae"}
#' (mean absolute error).
#' @param criterion The criterion used to select the tuning parameters. Available
#'  options are: \code{criterion = "cv"} which perform cross-validation to select
#'  the optimal tuning parameter, while \code{criterion = "aic"} and
#'  \code{criterion = "bic"} select the tuning parameter based on the AIC and BIC,
#'  respectively.
#' @param sigma2 The residual variance obtained by fitting the full OLS model without
#' any regularization or feature selection. It is used for the calculation of AIC and BIC.
#' To compute these metrics, the residual variance needs to be provided via the
#' `sigma2` parameter.
#' @param lambda.min.ratio When lambda values are automatically generated, the sequence
#' is determined by lambda.max and lambda.min ratio. The program generates nlambda values
#' on the log scale from lambda.max down to lambda.min. lambda.max is not user-specified but is
#' computed from the input standardized x and y. It is the smallest value for lambda such that
#' all the coefficients are zero. The default is lambda.min.ratio = 0.0001.
#' @param parallel The program use \code{"[cv.glmnet()]"} function which supports
#' parallel computing. To make it work, users must register parallel beforehand.
#' see `glmnet` package for details
#' @param betatypes Not used but added for consistency with the oracle_model().
#' @return Returns the following items:
#' \item{beta}{Shrunken regression estimates.}
#' \item{lambda}{The tuning parameter used for the estimation of parameters.}
#' \item{gamma}{Not applicable for ridge but returned for consistency with other methods like adaptive lasso.}
#' \item{x}{A standardized matrix of predictors used in fitting the linear model.}
#' \item{y}{A centered vector of the response variable used in fitting the linear model.}
#' @import glmnet
#' @export
ridgefit <- function(x, y,
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

  # assert that residual variance must be provided for AIC and BIC criteria
  if (criterion != "cv" && is.null(sigma2)) {
    stop(sprintf("!sigma2 must be provided for %s calculation.", toupper(criterion)))
  }

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # x must have column names
  if (is.null(colnames(x))) {
    stop("x must have column names.")
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

  # Generate sequence of lambdas when not provided. x and y assumed to be standardized
  if (is.null(lambda)) {
    lambda <- lambda_seq_ridge(x = x, y = y,
                               nlambda = nlambda,
                               epsilon = lambda.min.ratio
                              )
  }

  # override the default nlambda when lambdas are supplied by the user
  nlambda <- length(lambda)

  # regression coefficients without intercept for selected criterion. we assume
  # x and y are standardized especially if lambda was estimated by the program
  fit <- switch (criterion,
                  "cv" = ridge_cv(x = x, y = y,
                                  nfolds = nfolds,
                                  foldid = foldid,
                                  lambda = lambda,
                                  nlambda = nlambda,
                                  parallel = parallel,
                                  standardize = FALSE,
                                  type.measure =  type.measure
                                 ),

                  "aic" = ridge_aic_bic(x = x, y = y,
                                        lambda = lambda,
                                        nlambda = nlambda,
                                        standardize = FALSE,
                                        criterion = "aic",
                                        sigma2 = sigma2
                                       ),

                  "bic" = ridge_aic_bic(x = x, y = y,
                                        lambda = lambda,
                                        nlambda = nlambda,
                                        standardize = FALSE,
                                        criterion = "bic",
                                        sigma2 = sigma2
                                       )
  )
  # return regression estimates, x and y used to fit the model
  fit_parms <- list(beta = fit$beta,
                    lambda = fit$lambda,
                    gamma = NA,
                    x = x,
                    y = y
                   )

  # set class to be used by generic functions
  class(fit_parms) = "ridgefit"

  return(fit_parms)
}

#-----ridge with cross-validation-----------------------------------------------

#' Helper function that fits ridge with cross-validation tuning parameters
#' @inheritParams ridgefit
ridge_cv <- function(x, y,
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
  optimal_lambda <- glmnet::cv.glmnet(x = x, y = y,
                              foldid = foldid,
                              family = "gaussian",
                              lambda = lambda,
                              alpha= 0,
                              nlambda = nlambda,
                              type.measure = type.measure,
                              parallel = parallel,
                              standardize = standardize)$lambda.min

  fit <- glmnet::glmnet(x = x, y = y,
                lambda = optimal_lambda,
                family = "gaussian",
                standardize = standardize,
                alpha = 0
               )

  # regression estimates without the intercept
  out <- list(beta = setNames(as.numeric(fit$beta), colnames(x)),
              lambda = optimal_lambda)
  out
}

#-----ridge with AIC and BIC tuning parameters----------------------------------

#' Helper function that fits ridge with AIC and BIC tuning parameters
#' @inheritParams ridgefit
ridge_aic_bic <- function(x, y,
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
                alpha = 0,
                family = "gaussian"
               )

  # Matrix of regression coefficients without intercept
  betas <- fit$beta

  # linear predictor: crossprod in Matrix package supports sparse matrices
  xb <- as.matrix(Matrix::crossprod(t(x),betas))

  # residual for each value of lambda: each column is for a specific lambda
  #res <- as.matrix(-1*sweep(x = xb, MARGIN = 1, STATS = y, FUN = "-"))
  res <- matrix(y, nrow = nrow(xb), ncol = ncol(xb), byrow = FALSE) - xb

  # a vector of residual sum of squares for each value of lambda
  rss <- colSums(res ^ 2, na.rm = TRUE)

  # vector of degrees of freedom for each lambda: no intercept
  df <- df_ridge(x = x, lambda = lambda)

  # compute AIC and BIC for each value of lambda
  pn <- (rss / (n * sigma2))
  aic <- pn + (2 / n) * df
  bic <- pn + (log(n) / n) * df

  # Find lambda that yield smallest AIC and BIC
  lambda.aic <- lambda[which.min(aic)]
  lambda.bic <- lambda[which.min(bic)]

  # beta <- switch(criterion,
  #                "aic" = as.numeric(predict(fit,
  #                                           type = "coefficients",
  #                                           s = lambda.aic,
  #                                           exact = TRUE,
  #                                           alpha = 0,
  #                                           standardize = standardize,
  #                                           x = x,
  #                                           y = y))[-1],
  #                "bic" = as.numeric(predict(fit,
  #                                           type = "coefficients",
  #                                           s = lambda.bic,
  #                                           exact = TRUE,
  #                                           alpha = 0,
  #                                           standardize = standardize,
  #                                           x = x,
  #                                           y = y))[-1]
  #                )

  beta <- switch(criterion,
                 "aic" = as.numeric(glmnet(x = x, y = y,
                                           lambda = lambda.aic,
                                           standardize = standardize,
                                           alpha = 0,
                                           family = "gaussian")$beta),
                 "bic" = as.numeric(glmnet(x = x, y = y,
                                           lambda = lambda.bic,
                                           standardize = standardize,
                                           alpha = 0,
                                           family = "gaussian")$beta)
  )

  out <- list(beta = setNames(beta, colnames(x)),
              lambda = ifelse(criterion=="aic",lambda.aic,lambda.bic))
  return(out)
}

## see https://glmnet.stanford.edu/articles/glmnet.html
#' Function that calculates lambda sequence for the ridge
#' @param x a standardized matrix of predictors
#' @param y a centered vector of response
#' @param nlambda Number of lambdas
#' @param epsilon The ratio of smallest value of the generated lambda
#' sequence (say lambda.min) to lambda.max. See lambda.min.ratio
#' in \code{[cv.glmnet]} for more details
#' @export
lambda_seq_ridge <- function(x, y,nlambda = 100,epsilon =  0.0001){
  # calculate lambda-max that gives all estimates almost equal to 0
  # we use alpha = 0.001 instead of 0,
  # see https://glmnet.stanford.edu/articles/glmnet.html
  #lambda_max <- max(abs(t(x)%*%y))/(length(y)*0.001)
  #lambda_max <- max(abs(colSums(as.matrix(x) * y)), na.rm = T) / (length(y)*0.001)

  x <- Matrix::Matrix(as.matrix(x),sparse = TRUE)
  lambda_max <- max(abs(Matrix::crossprod(x, y)), na.rm = TRUE) / (length(y) * 0.001)

  # calculate lambda-min
  lambda_min <- epsilon*lambda_max

  # choose the grid of lambdas to be equidistant on the log-scale
  lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

  # return vector of lambdas
  return(lambda)
}

#' Extract coefficients from a ridgefit object
#' @method coef ridgefit
#' @param object Fitted \code{"ridgefit"} model object
#' @param ... not used at the moment
#' @return regression estimates without intercept
#' @export
coef.ridgefit <- function(object, ...) {
  object$beta
}

#' make predictions from a "ridgefit" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"ridgefit"} object.
#' @param object Fitted \code{"ridgefit"} model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a
#'@param ... not used at the moment
#' @method predict ridgefit
#' @export
predict.ridgefit <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  #newx = cbind(rep(1,dim(newx)[1L]),newx)
  newx %*% coef.ridgefit(object)
}


