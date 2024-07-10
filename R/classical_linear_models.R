#' @title Fit the Gaussian linear regression model with and without shrinkage.
#' @description Fits a linear regression model for a Gaussian response. The regression 
#' estimates can be subjected to shrinkage obtained from three different approaches: global
#'  shrinkage, parameter-wise shrinkage (PWS), and quadratic PWS (Breiman method). 
#' @param x A standardized matrix of predictors with dimensions nobs x nvars, where 
#' each row represents an observation vector. If not standardized, set standardize to
#' TRUE.
#' @param y A numeric, centered quantitative vector of responses. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y.
#' @param standardize Specifies whether to standardize the predictors matrix 
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. By default,  \code{standardize = FALSE}, indicating that no 
#'  standardization will be performed. This assumes that users have already standardized x 
#'  and centered y so that the intercept is zero.
#' @param foldid A vector of values between 1 and nfolds identifying what fold
#' each observation is in while conducting cross-validation. Default is NULL
#' and the program will generate it.
#' @param nfolds Number of folds for cross-validation. Default is 10. If k = n then k-fold is
#' equal to leave-one-out cross-validation.
#' @param shrinkage Specifies the method for shrinking parameter estimates.
#' Three options are available: Parameter-wise shrinkage (PWS), global shrinkage, and 
#' quadratic PWS (Breiman method). The default is `none`, and no shrinkage is conducted.
#' @param choice Specifies which method of estimating shrinkage factors to use: 
#' tenfold cross-validation or leave-one-out cross-validation. Default is 
#' tenfold cross-validation.
#' @param nonnegative Specify whether nonnegative shrinkage factors should be 
#' estimated. The default is FALSE, which means that negative shrinkage factors may be
#'  estimated. See the "lower.limits" argument for alternative instructions on estimating
#'   nonnegative shrinkage factors for the Breiman method.
#' @param lower.limits Vector of length nvar that specifies the lower limits for the shrinkage
#' factor when using Breiman's method. The default is -Inf, potentially allowing estimation of 
#' negative shrinkage factors. If \code{nonnegative = TRUE}, the program will set 
#' \code{lower.limits = rep(0, nvar)} for variables, leading to the estimation of nonnegative 
#' shrinkage factors. In short, setting  \code{lower.limits = 0} is identical to setting 
#' \code{nonnegative = TRUE}. This option is only applicable to the Breiman method.
#' @param upper.limits Vector of length nvar specifying upper limits for the shrinkage factors. 
#' The default is Inf, and unbounded shrinkage factors will be estimated. This 
#' option is only applicable to the Breiman method.
#' @param type.measure Loss to use for cross-validation. Currently, two options
#' are available. The default is `type.measure = "mse"`, which uses squared error. 
#' Another option is `type.measure = "mae"` (mean absolute error). This
#' option is only applicable to the Breiman method.
#' @param sigma2 Residual variance obtained by fitting the full model(model without variable
#'  selection). This value is utilized in the calculation of AIC and BIC, although it is 
#'  currently not in use.
#' @param betatypes Not used but added for consistency with the oracle_model().
#' @return: Returns a list with the following components:
#' \item{beta:}{Regression estimates. If shrinkage = "none", unshrunken OLS 
#' estimates are returned; otherwise, shrunken regression estimates are returned.}
#' \item{shrinkageFactors:}{Shrinkage factors for each predictor estimated 
#' using one of the specified shrinkage methods. If shrinkage = "none", the
#' shrinkage factors of 1s are returned.}
#' \item{nvar:}{The number of variables in the model.}
#' \item{x:}{A standardized matrix of predictors used in fitting the linear model.}
#' \item{y:}{A centered vector of the response variable used in fitting the 
#' linear model.}
#' \item{lambda,gamma:}{Parameters used in penalized methods, such as adaptive 
#' lasso. It is not applicable in the current context. However, the function 
#' returns the two parameters for consistency with the interface of penalized methods.}
#' @export 
linear_model <- function(x,
                         y,
                         foldid = NULL,
                         nfolds = 10,
                         shrinkage = c("none", "pwsf", "global", "breiman"),
                         choice = c("tenfold", "loocv"),
                         lower.limits = -Inf,
                         upper.limits = Inf,
                         type.measure = c("mse", "mae"),
                         nonnegative = FALSE,
                         standardize = FALSE,
                         sigma2 = NULL,
                         betatypes = NULL) { 
  
  # Match arguments
  shrinkage <-  match.arg(shrinkage)
  choice <-  match.arg(choice)
  type.measure <- match.arg(type.measure)
  
  # Initial estimate for Breiman's approach
  initial.estimates <-  "default"
  
  # Override lower limits when nonnegative is TRUE and the Breiman method is chosen.
  if (shrinkage == "breiman" && nonnegative) {
    lower.limits <- 0
  }
  
  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # Assert that x must have column names
  varx <- colnames(x)
  
  if (is.null(varx)){
      stop("x must have column names.")
  }
  
  
  # Sample size and number of variables
  dimx <- dim(x)
  n <- as.integer(dimx[1])
  nvar <- as.integer(dimx[2L])
  
  # Deal with foldid
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }
  
  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and column means
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)
    
    # center x, and scale x by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)
    
    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }
  
  # Fit linear model for standardized x and centered y
  beta <- .lm.fit(x = x, y = y)$coefficients

  # Default shrinkage factors for all variables in x
  shrinkageFactors <- setNames(rep(1, nvar), varx)
  
  # Shrink regression estimates based on shrinkage method and choice of tuning
  if (shrinkage!="none") {
    shrinkfit <- switch(shrinkage,
                        "pwsf" = switch(choice,
                                        "tenfold" = pwsf(x = x, y = y, 
                                                         varx = varx,
                                                         nfolds = nfolds,
                                                         foldid = foldid, 
                                                         loocv = FALSE, 
                                                         nonnegative = nonnegative,
                                                         standardize = FALSE
                                                         ),
                                        "loocv" =   pwsf(x = x, y = y,
                                                         varx = varx, 
                                                         nfolds = nfolds, 
                                                         foldid = foldid,
                                                         loocv = TRUE, 
                                                         nonnegative = nonnegative,
                                                         standardize = FALSE)
                        ),
                        "global" = switch(choice,
                                          "tenfold" = global(x = x, y = y,
                                                             varx = varx,
                                                             nfolds = nfolds, 
                                                             foldid = foldid,
                                                             nonnegative = nonnegative,
                                                             standardize = FALSE
                                                             ),
                                          "loocv" =   global(x = x, y = y, 
                                                             varx = varx, 
                                                             nfolds = n, # k = n is loocv
                                                             foldid = foldid,
                                                             nonnegative = nonnegative,
                                                             standardize = FALSE)
                        ),
                        "breiman" = switch(choice,
                                           "tenfold" = breiman(x = x, y = y,
                                                               varx = varx, 
                                                               nfolds = nfolds, 
                                                               foldid = foldid, 
                                                               initial.estimates = initial.estimates, 
                                                               type.measure = type.measure,
                                                               lower.limits = lower.limits,
                                                               standardize = FALSE
                                                               ),
                                           "loocv"   = breiman(x = x, y = y,
                                                               varx = varx, 
                                                               nfolds = n, 
                                                               foldid = foldid,
                                                               initial.estimates = initial.estimates,
                                                               type.measure = type.measure,
                                                               lower.limits = lower.limits,
                                                               standardize = FALSE
                                                               )
                        )
    )
    # regression coefficients without intercept
    beta <- shrinkfit$ShrunkenRegCoef
    
    # shrinkage factors without intercept
    shrinkageFactors <- shrinkfit$shrinkageFactors
  }
  
  fit <- list(
    beta = beta, # regression estimates without intercept: intercept assumed to be 0
    shrinkageFactors = shrinkageFactors,
    nvar = length(varx),
    lambda = NA,
    gamma = NA,
    x = x,
    y = y
  )

  class(fit) = "linearmodel"
  return(fit)
}


#' Extract coefficients from a "linearmodel" object
#'
#' @param object Fitted \code{"linearmodel"} model object
#' @param ... not used at the moment
#' @method coef linearmodel
#' @export
coef.linearmodel <- function(object, ...) {
  object$beta
}

#' @title Make predictions from a "linearmodel" object.
#'
#' @description
#' Make predictions from linear models of class \code{"linearmodel"}. It is 
#' similar to other predict methods in R.
#' @param object A fitted \code{"linearmodel"} model object.
#' @param newx A matrix of new values for \code{x} at which predictions are to 
#' be made. It must be a matrix.
#'@param ... not used at the moment
#' @method predict linearmodel
#' @export
predict.linearmodel <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  newx %*% coef.linearmodel(object)
}

