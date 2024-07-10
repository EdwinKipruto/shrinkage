#' @title Fit a linear regression model with best subset selection and backward
#' elimination method.
#' @description Conduct variable selection using best subset selection
#' (exhaustive search) and backward elimination. The final model can be chosen
#' using cross-validation, AIC, or BIC. The regression estimates
#' can be subjected to shrinkage, which is estimated using PWSF, global,
#' or the Breiman method. The function returns unshrunken or shrunken estimates.
#' @param x The standardized matrix of predictors.
#' @param y A numeric centered quantitative vector of response.
#' @param method Specifies which method to use for variable selection: either
#' exhaustive search or backward elimination method. The default is exhaustive
#' search via leaps and bound algorithm. See `leaps` package in R for details.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param foldid A vector of values between 1 and nfolds identifying which fold
#' each observation is in while conducting cross-validation. Used when criterion
#' is cross-validation (cv) or choice is tenfold. It can be missing when nfolds
#' is provided.
#' @param criterion Specifies which criterion to use when selecting the final
#' model. Default is cross-validation. Other available options are AIC and BIC.
#' @param shrinkage Specifies which shrinkage method of parameter estimates is
#' required. Three options are available. The default is `none`, and no shrinkage
#' is conducted. Parameter-wise shrinkage factors, global, or Breiman
#' method can be used to estimate shrinkage factors.
#' @param choice Specifies which method of estimating shrinkage factors to use:
#' tenfold cross-validation or leave-one-out cross-validation. Default is tenfold cv.
#' @param nonnegative Specify whether nonnegative shrinkage factors should be
#' estimated. The default is FALSE, implying thatnegative shrinkage factors may be
#'  estimated. See the `lower.limits` argument for instructions on estimating
#'  nonnegative shrinkage factors for the Breiman method.
#' @param lower.limits Vector of lower limits for each shrinkage factor.
#' The default is -Inf, and negative shrinkage factors will be estimated. If
#' \code{nonnegative = TRUE}, the program will set \code{lower.limits = 0}, and
#' nonnegative shrinkage factors will be estimated. In short, setting
#' \code{lower.limits = 0} is identical to setting \code{nonnegative = TRUE}.
#' This option is only applicable to the Breiman method.
#' @param upper.limits Vector of upper limits for the shrinkage factors.
#' The default is Inf, and unbounded shrinkage factors will be estimated. This
#' option is only applicable to the Breiman method
#' @param type.measure Loss to use for cross-validation. Currently two options,
#' are available. The default is type.measure="mse", which uses squared-error.
#' Other option include type.measure = "mae" (mean absolute error). This
#' option is only applicable to the Breiman method.
#' @param sigma2 Residual variance from the full model for calculation of AIC or
#' BIC. The user must supply if criterion is AIC or BIC.
#' @param standardize Specifies whether to standardize the predictors matrix
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. This standardization step can be useful to ensure
#'  that the predictors are on a comparable scale. By default,
#'  \code{standardize=FALSE}, indicating that no standardization will be performed.
#'  This assumes that users have already standardized their data.
#' @param lambda.min.ratio lambda can be provided if the user wants to specify
#' the lambda sequence, but typical usage is for the program to construct the
#' lambda sequence on its own. When automatically generated, the lambda sequence
#' is determined by lambda.max and lambda.min.ratio. The latter is the ratio of
#' smallest value of the generated lambda sequence (say lambda.min) to lambda.max.
#' The program generates nlambda values linear on the log scale from lambda.max
#' down to lambda.min. lambda.max is not user-specified but is computed from the
#' input x and y: it is the smallest value for lambda such that all the
#' coefficients are zero. Default is lambda.min.ratio = 1e-08.This option is only
#' applicable to the Breiman method.
#' @param betatypes Not used but added for consistency with the oracle_model().
#' @return Returns a list with the following components:
#' \item{beta:}{Regression estimates of standardized variables. If shrinkage is
#' not "none," shrunken estimates are returned.}
#' \item{nvar:}{The number of variables selected.}
#' \item{shrinkagefactors:}{Shrinkage factors for the selected variables.}
#' \item{x:}{Standardized design matrix of the selected variables.}
#' \item{y:}{Centered response variable.}
#' \item{lambda,gamma:}{Parameters used in penalized methods, such as adaptive
#' lasso. It is not applicable in the current context. However, the function
#' returns the two parameters for consistency with the interface of penalized methods.
#' }
#' @export
bestsubsetfit <- function(x, y,
                          method = c("exhaustive", "backward"),
                          foldid = NULL, nfolds = 10,
                          criterion = c("cv", "aic", "bic"),
                          shrinkage = c("none", "pwsf", "global", "breiman"),
                          choice = c("tenfold", "loocv"),
                          standardize = FALSE,
                          type.measure = c("mse", "mae"),
                          lower.limits = -Inf,
                          upper.limits = Inf,
                          nonnegative = FALSE,
                          lambda.min.ratio = 1e-08,
                          sigma2 = NULL,
                          betatypes = NULL){

  # Match arguments
  method <- match.arg(method)
  criterion <-  match.arg(criterion)
  shrinkage <-  match.arg(shrinkage)
  choice <-  match.arg(choice)

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # x must have column names
  xnames <- colnames(x)

  if (is.null(xnames)) {
    stop("x must have column names.")
  }

  # sample size and number of variables
  dimx <- dim(x)
  n <- as.integer(dimx[1L])
  nvar <- as.integer(dimx[2L])

  # set defaults
  if (shrinkage == "breiman" && nonnegative) {
    lower.limits <- 0
  }

  # assert that residual variance must be provided for AIC and BIC criteria
  if (criterion != "cv" && is.null(sigma2)) {
    stop(sprintf("!sigma2 must be provided for %s calculation.", toupper(criterion)))
  }

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of matrix x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)

    # center x and scale x by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }

  # best model of each size using exhaustive search or backward elimination
  model_fit <- leaps::regsubsets(x = x,
                                 y = y,
                                 method = method,
                                 nvmax = nvar
                                 )

  # Required for k-fold cross-validation
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # select the best model using CV, AIC or BIC
  if (criterion=="cv") {
    cvfit <- cv.regsubsets(x = x,
                           y = y,
                           nfolds = nfolds,
                           foldid = foldid,
                           method = method,
                           standardize = FALSE
                           )

    # we assume x is standardized and y centered so intercept = 0
    regcoef <- coef(model_fit, cvfit$bestmodel)[-1]
  } else {
    bfit <- regsubsetAICBIC(x = x, y = y, method = method, standardize = FALSE)
    regcoef <- switch (criterion,
                       "aic" = coef(model_fit, bfit$nvarAIC)[-1],
                       "bic" = coef(model_fit, bfit$nvarBIC)[-1]
    )
  }

  # names and number of variables selected
  varx <- names(regcoef)
  nx <- length(varx)

  # return betas like penalized methods: Eliminated variables padded with zeros
  index <- match(varx, xnames)
  beta <- setNames(rep(0, nvar), xnames)

  # Replace values at specified indices with regression coefficients
  beta[index] <- regcoef

  # shrinkage=="none" returns unshrunken estimates
  shrinkageFactors <- setNames(rep(0, nvar), xnames)
  shrinkageFactors[index] <- rep(1, nx)

  # shrunk regression estimates based on shrinkage method and choice of tuning
  if (shrinkage!="none") {
  shrinkfit <- switch(shrinkage,
               "pwsf" = switch(choice,
                               "tenfold" = pwsf(x = x,
                                                y = y,
                                                varx = varx,
                                                nfolds = nfolds,
                                                foldid = foldid,
                                                loocv = FALSE,
                                                nonnegative = nonnegative,
                                                standardize = FALSE
                                               ),
                               "loocv" =   pwsf(x = x,
                                                y = y,
                                                varx = varx,
                                                nfolds = nfolds,
                                                foldid = foldid,
                                                loocv = TRUE,
                                                nonnegative = nonnegative,
                                                standardize = FALSE)
               ),
               "global" = switch(choice,
                                 "tenfold" = global(x = x,
                                                    y = y,
                                                    varx = varx,
                                                    nfolds = nfolds,
                                                    foldid = foldid,
                                                    nonnegative = nonnegative,
                                                    standardize = FALSE
                                                   ),
                                 "loocv" =   global(x = x,
                                                    y = y,
                                                    varx = varx,
                                                    nfolds = n,
                                                    foldid = foldid,
                                                    nonnegative = nonnegative,
                                                    standardize = FALSE)
               ),
               "breiman" = switch(choice,
                                  "tenfold" = breiman(x = x,
                                                      y = y,
                                                      varx = varx,
                                                      nfolds = nfolds,
                                                      foldid = foldid,
                                                      initial.estimates = "default",
                                                      type.measure = type.measure,
                                                      lower.limits = lower.limits,
                                                      upper.limits = upper.limits,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      standardize = FALSE
                                                     ),
                                  "loocv"   = breiman(x = x, y = y,
                                                      varx = varx,
                                                      nfolds = n,
                                                      foldid = foldid,
                                                      initial.estimates = "default",
                                                      type.measure = "default",
                                                      lower.limits = lower.limits,
                                                      upper.limits = upper.limits,
                                                      lambda.min.ratio = lambda.min.ratio,
                                                      standardize = FALSE)
               )
  )
  # regression coefficients without intercept
  beta <- shrinkfit$ShrunkenRegCoef

  # shrinkage factors without intercept
  shrinkageFactors <- shrinkfit$shrinkageFactors
  }

  fit <- list(beta = beta, # Without intercept
              shrinkageFactors = shrinkageFactors,
              nvar = length(varx),
              x = x, y = y,
              lambda = NA,
              gamma = NA)

  class(fit) = "bestsubsetfit"
  return(fit)
}

#' @title Select the best subset using an information criterion.
#' @description Conduct variable selection via exhaustive search or backward elimination and
#' return the best model using an information criterion. The model with the smallest
#' AIC or BIC among models of different sizes is considered the best model.
#' @param x The standardized matrix of predictors.
#' @param y A numeric centered quantitative vector of response.
#' @param method Specifies which method to use for variable selection: either
#' exhaustive search or backward elimination method. The default is exhaustive
#' search via leaps and bound algorithm.
#' @param standardize Specifies whether to standardize the predictors matrix
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. This standardization step can be useful to ensure
#'  that the predictors are on a comparable scale. By default,
#'  \code{standardize=FALSE}, indicating that no standardization will be performed.
#'  This assumes that users have already standardized their data.
#' @return It returns a list with the following components:
#' \item{AIC:}{The AIC for each model size, where the first element corresponds
#' to the model with one variable, and so on.}
#' \item{BIC:}{The BIC for each model size, where the first element corresponds
#' to the model with one variable, and so on.}
#' \item{AIC.min:}{The minimum AIC}
#' \item{BIC.min:}{The mininum BIC}
#' \item{nvarAIC:}{The number of variables selected by smallest AIC}
#' \item{nvarBIC:}{The number of variables selected by smallest BIC}
#' @import stats
#' @export
regsubsetAICBIC <- function(x, y,
                            method = c("exhaustive", "backward"),
                            standardize = FALSE){

  # Match arguments
  method <- match.arg(method)

  # set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # sample size and number of variables
  dimx <- dim(x)
  n <- as.integer(dimx[1L])
  nvar <- as.integer(dimx[2L])

  # x must have column names
  xnames <- colnames(x)

  if (is.null(xnames)) {
    stop("x must have column names.")
  }

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = T)

    # center x by mean and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }

  # fit subset selection and return best model of each size
  bestmodel <- leaps::regsubsets(x = x,
                                 y = y,
                                 method = method,
                                 nvmax = nvar)

  # initialize aic and bic
  AICx <- BICx <- setNames(numeric(nvar), 1:nvar)
  logn <- log(n)
  for (j in 1:nvar) {
    # Each model has different number of variables (model size)
    coefx <- coef(bestmodel ,id=j)
    namex <- names(coefx)[-1]

    # fit linear model for x standardized and y centered: no intercept needed
    fitm <- .lm.fit(x = x[,namex, drop = FALSE], y = y)

    # df of the model. Number of independent predictors + scale parameter. we
    # do not include intercept since its 0 due to centering y and standardizing x
    dfx <- fitm$rank + 1

    # rss is identical to the deviance of the model
    #rss <- sum(fitm$residuals^2, na.rm = TRUE)

    # estimate aic and bic: we use the concept of AIC() in R. if our model
    # included an intercept the results would have been identical
    AICx[j] <- -2 * logl(res = fitm$residuals, w = NULL) + 2 * dfx
    BICx[j] <- -2 * logl(res = fitm$residuals, w = NULL) + logn * dfx
  }

  # position of minimum AIC and BIC
   pos_aic <- which.min(AICx)
   pos_bic <- which.min(BICx)

   fit <- list(AIC = AICx, # all aic for each model size
               BIC = BICx,  # all bic for each model size
               AIC.min = AICx[pos_aic], # minimum aic
               BIC.min = BICx[pos_bic], # min bic
               nvarAIC = pos_aic,       # position of model that returned min aic
               nvarBIC = pos_bic        # position of model that returned min bic
              )
  return(fit)
}

## loglikehood formula for lm model. important if we decide to use lm.fit
# that does not include intercept
## source code: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/logLik.R
## aic = -2logL + 2*k where k is total number of parameters including
## intercept (if estimated) and residual variance. so if we use lm.fit
## intercept is not estimated because the data is assumed to be centered
## the total number of parameters will be less by one. compare results
## with AIC() instead of extractAIC()
## logl of lm and lm.fit are identical
## from the formula below: aic = -2*logl(model$residuals) + 2*k

#' estimates loglikehood of a Gaussian model
#'
#' @param res residuals from a fitted model. Should not contain missing values
#' @param w an optional vector of weights used in the fitting process. Should
#' be NULL or a numeric vector. see `[lm()]` for more details
logl <- function(res, w = NULL){
  N <- length(res)
  if (is.null(w)) {
    w <- rep.int(1, N)
  }
0.5* (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w*res^2))))
}

#' @title Cross-validation for selecting best model in subset selection
#' @description Select the best model using K-fold cross-validation. The mean squared error
#' (MSE) is utilized as the loss function. The model with the smallest average
#' MSE is considered the best model.
#' @param x A matrix of predictors with column names
#' @param y A vector of quantitative response variable.
#' @param nfolds The number of cross-validation folds.
#' @param foldid An optional vector identifying the fold each observation belongs to.
#' Values should range from 1 to k. If supplied, k can be missing.
#' @param method The method for variable selection, either "exhaustive" or "backward".
#' @param standardize Specifies whether to standardize the predictors matrix `x`
#' to have mean 0 and unit variance. If set to TRUE, `x` will be standardized.
#' Additionally, the response variable `y` will be centered by subtracting its
#' mean. This standardization step can be useful to ensure that the predictors
#' are on a comparable scale. By default, `standardize` is set to FALSE,
#' indicating that no standardization will be performed.
#' @references James, G., Witten, D., Hastie, T., & Tibshirani, R. (2013).
#' An Introduction to Statistical Learning (Vol. 112, p. 18). New York: Springer.
#'
#' @importFrom leaps regsubsets
#' @export
cv.regsubsets <- function(x, y,
                          nfolds = 10,
                          foldid = NULL,
                          method = c("exhaustive", "backward"),
                          standardize = F){

  # assert that x must have column names
  if (is.null(colnames(x))){
    stop("x must have column names")
  }

  # Match arguments
  method <- match.arg(method)

  # standardize x and center y if necessary
  if (standardize) {
    # sd and mean of x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = T)

    # center x by mean and scale by sd
    x <- sweep(x, 2L, colx_means, "-",  check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }

  # restrict number of k-folds to 3
  if (nfolds < 3) {
    stop("The number of folds must not be less than 3")
  }

  # sample size and number of variables
  dimx <- dim(x)
  n <- as.integer(dimx[1L])
  nvar <- as.integer(dimx[2L])

  # partitions the group into approx equal sizes
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # important if nfolds = n which is equivalent to leave-one-out cross-validation
  nfolds <- length(unique(foldid))

  # perform cv for each model size
  cv_errors <- matrix (NA, nfolds, nvar, dimnames = list(NULL, paste(1:nvar)))

  for (j in 1:nfolds) {
    # conduct exhaustive search or backward selection on the training set.
    best.fit <- leaps::regsubsets(x = x[foldid!=j,,drop=FALSE],
                                  y = y[foldid!=j],
                                  nvmax=nvar,
                                  method = method
                                 )

    # calculate MSE on the hold out data for each model size starting from
    # a model with one variable to a model with nvar variables
    for (i in 1:nvar) {
      pred <- predict(best.fit, x[foldid == j, , drop = FALSE], id = i)
      cv_errors[j,i]<- mean((y[foldid == j] - pred) ^ 2, na.rm = TRUE)
    }
  }
  # calculate the mean cross validation errors for each model size with their
  # corresponding standard deviation
  mean_cv_errors <- setNames(colMeans(cv_errors, na.rm = TRUE), seq(1, nvar))
  sd_cv_errors <- setNames(apply(cv_errors, 2, sd), seq(1, nvar))

  # standard error of mean_cv_errors = sd/sqrt(length(cv_errors))
  nn <- apply(cv_errors, 2, length)
  cvm_se <- setNames((sd_cv_errors / sqrt(nn)), seq(1, nvar))

  #upper curve = cvm + cvm_se and lower curve
  cvup <- mean_cv_errors + cvm_se
  cvlo <- mean_cv_errors - cvm_se

  # minimum cvm
  cvm_min <- min(mean_cv_errors, na.rm = TRUE)

  # The size of the best model
  bestmodel <- which.min(as.numeric(mean_cv_errors))

  # Return
  fit <- list(cv_errors = cv_errors,
              cvm = mean_cv_errors,
              cvsd = sd_cv_errors,
              cvm.se = cvm_se,
              cvup = cvup,
              cvlo = cvlo,
              cvm.min = cvm_min,
              bestmodel = bestmodel
              )
  class(fit) = "cv.regsubsets"
  return(fit)
}

#' @title Predict Method for regsubsets
#' @description
#' Make predictions from a "regsubsets" object. Similar to other predict methods.
#' @param object Fitted "regsubsets" model.
#' @param newx Matrix of new values for 'x' at which predictions are to be made. Must be a matrix.
#' @param id The model size at which predictions are required. Note that different model sizes are returned by regsubsets.
#' @param ... not used at the moment
#' @method predict regsubsets
#' @export
predict.regsubsets <- function (object,newx = NULL,id,...){

  # Check column names
  if (!is.null(newx) && is.null(colnames(newx))) {
    stop("newx must have names")
  }

  if (is.null(newx)) {
    stop("newx is NULL and must be provided")
  }

  # extract regression coefficients without intercept
  coefi <- coef(object,id = id)[-1]
  xvars <- names(coefi)

  newx[,xvars,drop = FALSE]%*%coefi
}


#' @title Extract coefficients from a bestsubsetfit object.
#' @description
#' It returns shrunken regression coefficients if `shrinkage = TRUE` in the `bestsubsetfit` object
#' @param object Fitted \code{"alassofit"} model object
#' @method coef bestsubsetfit
#' @param ... functionless
#' @export
coef.bestsubsetfit = function(object, ...) {
  object$beta
}

#' @title make predictions from a "bestsubsetfit" object.
#' @description
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"bestsubsetfit"} object.
#' @param object Fitted \code{"bestsubsetfit"} model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict bestsubsetfit
#' @export
predict.bestsubsetfit <- function(object, newx = NULL, ...) {

  if (is.null(newx)) {
    newx <- object$x
  }

  newx %*% coef.bestsubsetfit(object)
}
