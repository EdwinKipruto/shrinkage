#' @title Fit a Gaussian linear model with nonnegative garrote regularization
#'
#' @description Fits the nonnegative garrote regression for continuous outcomes proposed by
#' Breiman (1995). Currently, it does not support other types of outcomes such as categorical or
#' survival.
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
#' standardized to have unit variance. This is important if standardize.response was
#' set to TRUE in `cv.garrote()`. This option divides \code{`y`} by its standard deviation.
#' @param lambda  Regularization parameter that controls the strength of the penalty on shrinkage
#' factors.
#' @param lambda.initial The tuning parameter for the initial estimates using
#' lasso or ridge regularization. If \code{initial.estimates = "default"},
#' `lambda.initial` can be NULL since OLS estimates do not require a tuning
#' parameter. However, if \code{initial.estimates = "ridge"} or
#' \code{initial.estimates = "lasso"}, the tuning parameters must be supplied.
#' @param beta.initial A vector of initial estimates obtained from standardized
#' predictor matrix x. The default value is \code{beta.initial = NULL}, and the
#' program will estimate them based on the selected `initial.estimate` and
#' `lambda.initial`.
#' @param initial.estimates Specifies the type of initial estimates to use for
#' constructing nonnegative garrote weights. The \code{initial.estimates = "default"}
#' option uses \code{"OLS"} estimates. Other options include
#' \code{initial.estimates ="ridge"} or \code{initial.estimates ="lasso"}
#' estimates, which are calculated using the user-supplied \code{`lambda.initial`.}
#' @param alpha The mixing parameter where alpha = 1 represents the usual nonnegative
#' garrote, while alpha = 0 represents the garrote with a ridge-like penalty.
#' @param lower.limits Vector of lower limits for each shrinkage factor.
#' The default is 0, which will be replicated to have a length equal to the number of predictors, and
#' nonnegative shrinkage factors will be estimated.
#' @param upper.limits Vector of upper limits for the shrinkage factors. The default is Inf, and
#' unbounded shrinkage factors will be estimated.
#' @return Returns the following items:
#'  \item{beta}{Shrunken regression estimates without an intercept}
#'  \item{shrinkageFactors}{Shrinkage factors for all regression estimates}
#'  \item{lambda}{Lambda value used for tuning the model}
#'  \item{x}{Matrix \code{x} used to fit the model}
#'  \item{y}{Values of \code{y} used to fit the model}
#' @importFrom Matrix Matrix crossprod
#' @export
garrote <- function(x, y,
                    lambda,
                    lambda.initial = NULL,
                    beta.initial = NULL,
                    alpha = 1,
                    standardize = FALSE,
                    standardize.response = FALSE,
                    lower.limits = 0,
                    upper.limits = Inf,
                    initial.estimates = c("default", "ridge", "lasso")){

  # Match arguments
  initial.estimates <- match.arg(initial.estimates)

  # Assert that initial estimates must be provided if lambda.initial is missing
  if (initial.estimates != "default" && is.null(lambda.initial)) {
    stop("Provide a value for lambda.initial or set initial.estimates = 'default'")
  }

  # Check the length of lambda.initial
  nt <- length(lambda.initial)

  # If lambda.initial is not NULL and its length is not 1, raise an error
  if (!is.null(lambda.initial) && nt != 1) {
    stop(sprintf("The length of lambda.initial must be 1, not %d", nt))
  }

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # sample size and number of variables
  np <- dim(x)
  n <- as.integer(np[1L])
  nvar <- as.integer(np[2L])

  # Assign column names to variable xnames
  xnames <- colnames(x)

  # Check if xnames is NULL
  if (is.null(xnames)) {
    # If xnames is NULL, create default column names "x1", "x2", ..., "xnvar"
    xnames <- sprintf("x%d", 1:nvar)

    # Assign the generated column names to the columns of x
    colnames(x) <- xnames
  }

  # Check the length of lambda
  if (length(lambda) != 1) {
    # If length of lambda is not 1, raise an error
    stop("Only one value of lambda is currently supported")
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
    y <- y - mean(y, na.rm = TRUE)
  }

  # Standardize the response variable to have sd of 1
  if(standardize.response){
    y <- y / sd_x(y)
  }

  # Initial estimates for standardized x and centered/standardized y
  if (is.null(beta.initial)) {

    beta.initial = switch(initial.estimates,
                          "default" = stats::.lm.fit(x = x, y = y)$coefficients,
                          "ridge" = as.numeric(glmnet(x = x, y = y,
                                                      lambda = lambda.initial,
                                                      alpha = 0,
                                                      intercept = TRUE,
                                                      family = "gaussian",
                                                      standardize = FALSE)$beta),
                          "lasso" = as.numeric(glmnet(x = x, y = y,
                                                      lambda = lambda.initial,
                                                      alpha = 1,
                                                      intercept = TRUE,
                                                      family = "gaussian",
                                                      standardize = FALSE)$beta))
  }

  # Create a sparse matrix matb with diagonal elements from beta.initial
  matb <- Matrix::Matrix(diag(beta.initial,nrow = length(beta.initial)), sparse = TRUE)

  # Calculate the product x.tilde = x * beta.initial using cross product
  x.tilde <- Matrix::crossprod(Matrix(Matrix::t(x),sparse = TRUE), matb)

  # Estimate shrinkage factors
  shrinkage_factors <- as.numeric(glmnet::glmnet(x = x.tilde,
                                    y = y,
                                    lambda = lambda,
                                    standardize = FALSE,
                                    intercept = FALSE,
                                    lower.limits = lower.limits, # Important: c>=0
                                    upper.limits = upper.limits,
                                    family = "gaussian",
                                    weights = NULL,
                                    offset = NULL,
                                    alpha = alpha)$beta
                     )

  # Garrote standardized estimates: beta.initial*shrinkage_factors
  betas <- setNames(beta.initial * shrinkage_factors, xnames)
  names(shrinkage_factors) <- xnames

  # Return a list of parameters of interest
  fit <- list(
    beta = betas,    # no intercept
    shrinkageFactors = shrinkage_factors,
    lambda = lambda,
    x = x,
    y = y
  )

  class(fit)= "garrote"
  return(fit)
}


#' @title Extract coefficients from a garrote object
#'
#' @description
#' Extract the coefficients from a fitted garrote model. The regression coefficients
#' are provided without the intercept.
#'
#' @param object Fitted \code{"garrote"} model
#' @method coef garrote
#' @param ... Not used at the moment
#'
#' @return Returns the regression coefficients
#' @export
coef.garrote <- function(object, ...) {
  object$beta
}

#' @title make predictions from a "garrote" object.
#' @description
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"garrote"} object.
#' @param object Fitted \code{"garrote"} model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict garrote
#' @export
predict.garrote <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  newx %*% coef.garrote(object)
}

#' @title Cross-validation for nonnegative garrote
#'
#' @description
#' Perform k-fold cross-validation for nonnegative garrote and return the optimal value of
#' lambda. The estimated value of lambda depends on the initial estimate used.
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
#' standardized to have unit variance. This is important if standardize.response was
#' set to TRUE in `cv.garrote()`. This option divides \code{`y`} by its standard deviation.
#' @param lambda.initial An optional tuning parameter for the initial estimates
#' using lasso or ridge regularization. If a vector of lambda is provided,the
#' optimal will be determined by CV
#' @param beta.initial A vector of initial estimates obtained from standardized
#' predictor matrix x. The default value is \code{beta.initial = NULL}, and the
#' program will estimate them based on the selected `initial.estimate` and
#' `lambda.initial`.
#' @param initial.estimates Specifies the type of initial estimates to use for
#' constructing adaptive lasso weights. The \code{initial.estimates = "default"}
#' option uses \code{"OLS"} estimates. Other options include
#' \code{initial.estimates ="ridge"} or \code{initial.estimates ="lasso"}
#' estimates, which are calculated using the user-supplied \code{`lambda.initial`.}
#' @param nfolds The number of folds used when cross-validating lambda as well
#' as the lasso and ridge initial estimates. Default is 10.
#' @param foldid	an optional vector of values between 1 and nfold identifying
#' what fold each observation is in. If supplied, nfolds can be missing.
#' @param lambda Optional user-supplied sequence of lambdas. The default is NULL, and
#' the program chooses its own sequence. If lambda is supplied, its length must be at least two.
#' @param nlambda An integer specifying the number of values of lambda
#' for both nonnegative garrote and initial estimators (ridge and lasso). The default is NULL,
#' which sets `nlambda` to 100.
#' @param type.measure Loss function to use for cross-validation for both nonnegative garrote and
#' initial estimates (ridge or lasso) when not provided. Currently, two options are available. The
#' default option is \code{type.measure = "mse"}, which corresponds to squared error. Another
#' option is \code{type.measure = "mae"}, which corresponds to mean absolute error.
#' @param lambda.min.ratio Lambda can be provided if the user wants to specify
#' the lambda sequence, but the typical usage is for the program to construct the
#' lambda sequence on its own. When automatically generated, the lambda sequence
#' is determined by the ratio of lambda.max and lambda.min. The program generates nlambda values
#' on the log scale from lambda.max down to lambda.min. Lambda.max is not user-specified but is
#' computed from the input standardized x and y. It is the smallest value for lambda in which all
#' the estimated coefficients are zero. The default is `lambda.min.ratio = 0.0001`.
#' @param alpha The mixing parameter where alpha = 1 represents the usual nonnegative
#' garrote, while alpha = 0 represents the garrote with a ridge-like penalty.
#' @param lower.limits Vector of lower limits for each shrinkage factor.
#' The default is 0, which will be replicated to have a length equal to the number of predictors, and
#' nonnegative shrinkage factors will be estimated.
#' @param upper.limits Vector of upper limits for the shrinkage factors. The default is Inf, and
#' unbounded shrinkage factors will be estimated.
#' @return  Returns \code{cv.glmnet} objects. Please refer to the return
#' section of the mentioned function in the \code{glmnet} package. Additionally,
#' it returns the lambda (lambda.initial) for the optimal lasso or ridge initial
#' estimates.
#' @export
cv.garrote <- function(x, y,
                       nfolds = 10,
                       foldid = NULL,
                       lambda = NULL,
                       alpha = 1,
                       nlambda = 50,
                       lower.limits = 0,
                       upper.limits = Inf,
                       standardize = FALSE,
                       standardize.response = FALSE,
                       lambda.initial= NULL,
                       beta.initial = NULL,
                       initial.estimates = c("default", "ridge", "lasso"),
                       type.measure = c("mse", "mae"),
                       lambda.min.ratio =1e-04) {

  # match arguments
  initial.estimates <- match.arg(initial.estimates)
  type.measure <- match.arg(type.measure)

  # Ensure that the number of lambdas provided is 2 or more
  if (!is.null(lambda)) {
  # Count the number of lambdas
  num_lambdas <- length(lambda)

  # Check if the number of lambdas is less than 2
  if (num_lambdas < 2) {
    # If fewer than 2 lambdas, raise an error
    stop(sprintf("The number of lambdas (%d) must be 2 or more", num_lambdas))
  }
  }

  # set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # sample size and number of variables
  np <- dim(x)
  n <- as.integer(np[1])
  nvar <- as.integer(np[2])

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = T)

    # center x by mean and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y
    y <- y - mean(y, na.rm = TRUE)
  }

  # Standardize the response variable to have sd of 1
  if(standardize.response){
    y <- y / sd_x(y)
  }

  # Divide the data into approx equal groups for cross-validation
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # One-dimensional cross-validation. Initial estimates are not tuned; they are estimated
  # and then held constant across all lambda values. Assumed to be on the standardized scale.

  lambda.init <- lambda.initial
  beta.init <- beta.initial

  if(is.null(beta.init)){

    fitm <- switch(initial.estimates,
                         "default" = stats::.lm.fit(x = x, y = y),
                         "ridge" = glmnetCoef(x = x, y = y,
                                              foldid = foldid,
                                              lambda = lambda.initial,
                                              nlambda = nlambda,
                                              type.measure = type.measure,
                                              alpha = 0,
                                              standardize = FALSE
                                              ),
                         "lasso" = glmnetCoef(x = x, y = y,
                                              foldid = foldid,
                                              lambda = lambda.initial,
                                              nlambda = nlambda,
                                              type.measure = type.measure,
                                              alpha = 1,
                                              standardize = FALSE))
    # initial coefficients for predictors without intercept
    beta.init<- fitm$coefficients
    lambda.init <- ifelse(!is.null(fitm$lambda.initial), fitm$lambda.initial, 0)
  }

  # Generate sequence of lambdas when not provided
  if (is.null(lambda)) {
    lambda <- lambda_seq_garrote(x = x, y = y,
                                 beta.initial = beta.init,
                                 nlambda = nlambda,
                                 epsilon =  lambda.min.ratio,
                                 alpha = alpha
                                 )
  }

  # New values of x after multiplication with initial estimates. Utilize the
  # sparsity of resulting matrices for computational efficiency
  matb <- Matrix::Matrix(diag(beta.init,nrow = length(beta.init)), sparse = TRUE)

  x.tilde <- Matrix::crossprod(Matrix(Matrix::t(x),sparse = TRUE),matb)

  # Obtain the lambda that produces optimal shrinkage factors using cv.glmnet().
  cvfit <- glmnet::cv.glmnet(x = x.tilde,
                             y = y,
                             lambda = lambda,
                             standardize = FALSE,
                             intercept = FALSE,
                             foldid =  foldid, # important to avoid setting local seed
                             type.measure = type.measure,
                             alpha = alpha,
                             nlambda = nlambda,
                             lower.limits = lower.limits, # Important: c>=0
                             upper.limits  = upper.limits,
                             lambda.min.ratio = lambda.min.ratio,
                             family = "gaussian"
                             )

  # add lambda for initial estimates to be used by garrote()
  cvfit$lambda.initial <- lambda.init

  return(cvfit)
}

#' @title Fit a linear model with nonnegative garrote regularization
#'
#' @description Fits the nonnegative garrote model proposed by Breiman (1995) using the optimal
#' tuning parameter determined through cross-validation. It also supports the selection
#' of tuning parameters based on the Akaike Information Criterion (AIC) or the Bayesian Information
#' Criterion (BIC). The value of lambda that yields the smallest AIC or BIC is chosen as the best
#' tuning parameter.
#' @param x A standardized matrix of predictors of dimension nobs x nvars, where
#' each row is an observation vector.If not standardized, set standardize to TRUE.
#' @param y A numeric centered quantitative vector of response. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. See `standardize.response`
#' for more details.
#' @param standardize Specifies whether to standardize the predictors matrix
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#' standardized. Additionally, the response variable \code{`y`} will be centered
#' by subtracting its mean. This standardization step can be useful to ensure
#' that the predictors are on a same scale. By default,
#' \code{standardize = FALSE}, indicating that no standardization will be performed.
#' This assumes that users have already standardized x and  y so that the intercept
#' is zero.
#' @param standardize.response Specifies whether the response variable y should be
#' standardized to have unit variance. This option divides \code{`y`} by its standard deviation.
#' @param criterion The criterion used to select tuning parameters for the
#' nonnegative garrote. Available options include: \code{criterion = "cv"}
#' which perform cross-validation to select the optimal tuning parameter, while
#' \code{criterion = "aic"} and \code{criterion = "bic"} select the tuning parameter
#'  based on the AIC and BIC, respectively.
#' @param lambda Optional user-supplied sequence of lambdas. The default is NULL, and
#' the program chooses its own sequence. If lambda is supplied, its length must be at least two.
#' @param nlambda An integer specifying the number of values of lambda
#' for both nonnegative garrote and initial estimators (ridge or lasso). The default is NULL,
#' which sets `nlambda` to 100.
#' @param initial.estimates Specifies the type of initial estimates to use for
#' constructing nonnegative garrote weights. The \code{initial.estimates = "default"}
#' option uses \code{"OLS"} estimates. Other options include
#' \code{initial.estimates ="ridge"} or \code{initial.estimates ="lasso"}
#' estimates, which are calculated using the user-supplied \code{`lambda.initial`.}
#' @param foldid an optional vector of values between 1 and nfolds identifying
#' what fold each observation is in. If supplied, nfolds can be missing.
#' @param nfolds Number of folds for cross-validation. Default is 10. Smallest value allowable is
#' nfolds = 3.
#' @param type.measure Loss function to use for cross-validation for both nonnegative garrote and
#' initial estimates (ridge or lasso) when not provided. Currently, two options are available. The
#' default option is \code{type.measure = "mse"}, which corresponds to squared error. Another
#' option is \code{type.measure = "mae"}, which corresponds to mean absolute error.
#' @param lambda.min.ratio Lambda can be provided if the user wants to specify
#' the lambda sequence, but the typical usage is for the program to construct the
#' lambda sequence on its own. When automatically generated, the lambda sequence
#' is determined by the ratio of lambda.max and lambda.min. The program generates nlambda values
#' on the log scale from lambda.max down to lambda.min. Lambda.max is not user-specified but is
#' computed from the input standardized x and y. It is the smallest value for lambda in which all
#' the estimated coefficients are zero. The default is `lambda.min.ratio = 0.0001`.
#' @param sigma2 The residual variance obtained by fitting the full OLS model without
#' any regularization or feature selection. It is used for the calculation of AIC and BIC.
#' To compute these metrics, the residual variance needs to be provided via the
#' `sigma2` parameter.
#' @param alpha The mixing parameter where alpha = 1 represents the usual nonnegative
#' garrote, while alpha = 0 represents the garrote with a ridge-like penalty.
#' @param lower.limits Vector of lower limits for each shrinkage factor.
#' The default is 0, which will be replicated to have a length equal to the number of predictors, and
#' nonnegative shrinkage factors will be estimated.
#' @param upper.limits Vector of upper limits for the shrinkage factors. The default is Inf, and
#' unbounded shrinkage factors will be estimated.
#' @param betatypes Not used but added for consistency with the oracle_model().
#' @param use_alasso_lambda Specify whether the adaptive lasso sequence of lambdas should be used.
#' This depends on the value of `gamma` supplied. See `lambda_seq_alasso()`. Only calculated when
#' lambda = NULL
#' @param gamma adaptive lasso parameter used to calculate lambdas when use_alasso_lambda is set
#' to TRUE. Default is gamma = 1
#' @param ignore.gamma see lambda_seq_alasso() for details
#' @return returns the following items:
#' \item{beta}{Shrunken regression estimates.}
#' \item{shrinkageFactors}{Shrinkage factors for each predictor.}
#' \item{lambda}{The tuning parameter used for estimation of parameters.}
#' \item{gamma}{Not applicable for garrote but returned for consistency with other methods like adaptive lasso.}
#' \item{x}{A standardized matrix of predictors used in fitting the linear model.}
#' \item{y}{A centered vector of the response variable used in fitting the linear model.}
#' @export
nngfit <- function(x, y,
                   standardize = FALSE,
                   standardize.response = FALSE,
                   criterion = c("cv", "aic", "bic"),
                   lambda = NULL,
                   nlambda = 100,
                   initial.estimates = c("default", "ridge", "lasso"),
                   nfolds = 10,
                   foldid = NULL,
                   type.measure = c("mse", "mae"),
                   lambda.min.ratio = 1e-04,
                   alpha = 1,
                   sigma2 = NULL,
                   lower.limits = 0,
                   upper.limits = Inf,
                   use_alasso_lambda = FALSE,
                   ignore.gamma = FALSE,
                   gamma = 1,
                   betatypes = NULL){

  # Match arguments
  initial.estimates <- match.arg(initial.estimates)
  type.measure <- match.arg(type.measure)
  criterion <- match.arg(criterion)

  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # assert that residual variance must be provided for AIC and BIC criteria
  if (criterion != "cv" && is.null(sigma2)) {
    stop(sprintf("!sigma2 must be provided for %s calculation.", toupper(criterion)))
  }

  # sample size and number of variables
  np <- dim(x)
  n <- as.integer(np[1])
  nvar <- as.integer(np[2])

  # standardize x and center y if necessary: important if intercept is not needed
  if (standardize) {
    # standard deviation and mean of x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)

    # center x and scale x by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)

    # center the response y: DIVIDED BY SD TO RETURN CORRECT LAMBDAS
    y <- y - mean(y, na.rm = TRUE)
  }

  # standardize y
   if (standardize.response) {
     y <- y / sd_x(y)
   }

  # a vector of values between 1 and nfold used to split data in cv
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }

  # Obtain initial estimates: OLS, ridge, or lasso.
  fitn <- switch(initial.estimates,
                       "default" = stats::.lm.fit(x = x, y = y),
                       "ridge" = glmnetCoef(x = x, y = y,
                                            foldid = foldid,
                                            nlambda = nlambda,
                                            type.measure = type.measure,
                                            alpha = 0,
                                            standardize = FALSE,
                                            lambda = NULL
                                            ),
                       "lasso" = glmnetCoef(x = x, y = y,
                                            foldid = foldid,
                                            nlambda = nlambda,
                                            type.measure = type.measure,
                                            alpha = 1,
                                            standardize = FALSE,
                                            lambda = NULL)
                 )

  # Initial estimates and lambda that was used to obtain the estimates.
  beta.init <- fitn$coefficients

  lambda.init<- ifelse(!is.null(fitn$lambda.initial),fitn$lambda.initial, 0)

  # sequence of lambdas when not provided
  if (is.null(lambda)) {
    # We can use lambdas from adaptive lasso for comparison if needed
    if (use_alasso_lambda){
      lambda <- lambda_seq_alasso(x = x, y = y,
                                  beta.initial = beta.init,
                                  nlambda = nlambda,
                                  epsilon =  lambda.min.ratio,
                                  gamma = gamma,
                                  ignore.gamma = ignore.gamma)
    } else{
    lambda <- lambda_seq_garrote(x = x, y = y, # Y SHOULD BE DIVIDED BY SD SINCE ITS CENTERED, NO PROBLEM FOR OUR SIMULATION
                                 beta.initial = beta.init,
                                 nlambda = nlambda,
                                 epsilon =  lambda.min.ratio,
                                 alpha = alpha)
    }
  }

  # The number of lambda values can change from the default when the user supplies lambdas.
  nlambda <- length(lambda)

  # fit models using different criteria
  fit <- switch (criterion,
                 "cv" = nng_cv(x = x, y = y,
                               nfolds = nfolds,
                               foldid = foldid,
                               lambda = lambda,
                               type.measure = type.measure,
                               nlambda = nlambda,
                               beta.initial = beta.init,
                               lambda.initial = lambda.init,
                               initial.estimates = initial.estimates,
                               lower.limits = lower.limits,
                               upper.limits = upper.limits,
                               alpha = alpha,
                               standardize = FALSE,
                               standardize.response = FALSE,
                               lambda.min.ratio =lambda.min.ratio
                               ),
                 "aic" = nng_aic_bic(x = x, y = y,
                                     beta.initial = beta.init,
                                     initial.estimates = initial.estimates,
                                     lambda.initial = lambda.init,
                                     lower.limits = lower.limits,
                                     upper.limits = upper.limits,
                                     criterion = "aic",
                                     lambda = lambda,
                                     alpha = alpha,
                                     sigma2 = sigma2
                                     ),
                 "bic" = nng_aic_bic(x = x, y = y,
                                     beta.initial = beta.init,
                                     initial.estimates = initial.estimates,
                                     lambda.initial = lambda.init,
                                     lower.limits = lower.limits,
                                     upper.limits = upper.limits,
                                     criterion = "bic",
                                     alpha = alpha,
                                     lambda = lambda,
                                     sigma2 = sigma2
                                     )
  )

  # Return betas, lambda, x, y, and gamma. Gamma is not applicable in nng, but
  # it is used in alasso and relaxed lasso. To maintain uniformity, we return NA.
  fitx <- list(
      beta = fit$beta,
      shrinkageFactors = fit$shrinkageFactors,
      lambda = fit$lambda,
      gamma = NA,
      x = x,
      y = y
    )

  # set class for generic functions
  class(fitx) = "nngfit"
  return(fitx)
}

#' Helper function for `nngfit()`
#'
#' Helper function that fits nonnegative garrote using optimal tuning parameters
#' from cross-validation
#' @inheritParams nngfit
#' @param lambda.initial The tuning parameter for the initial estimates using
#' lasso or ridge regularization. If \code{initial.estimates = "default"},
#' `lambda.initial` can be NULL since OLS estimates do not require a tuning
#' parameter. However, if \code{initial.estimates = "ridge"} or
#' \code{initial.estimates = "lasso"}, the tuning parameters must be supplied.
#' @param beta.initial A vector of initial estimates obtained from standardized
#' predictor matrix x. The default value is \code{beta.initial = NULL}, and the
#' program will estimate them based on the selected `initial.estimate` and
#' `lambda.initial`.
nng_cv <- function(x, y, nfolds,
                   foldid, lambda,
                   type.measure,
                   nlambda,
                   beta.initial,
                   lambda.initial,
                   initial.estimates,
                   lower.limits,
                   upper.limits,
                   standardize,
                   standardize.response,
                   lambda.min.ratio,
                   alpha){

  # Estimate optimal lambda via cross-validation
  optimal_lambda <- cv.garrote(x = x, y = y,
                               nfolds = nfolds,
                               foldid = foldid,
                               lambda = lambda,
                               beta.initial = beta.initial,
                               type.measure = type.measure,
                               nlambda = nlambda,
                               lambda.initial = lambda.initial,
                               alpha = alpha,
                               initial.estimates = initial.estimates,
                               lower.limits = lower.limits,
                               upper.limits = upper.limits,
                               standardize = standardize,
                               standardize.response = standardize.response,
                               lambda.min.ratio = lambda.min.ratio)$lambda.min

  # Fit nonnegative garrote model using optimal lambda
  fit <- garrote(x = x, y = y,
                 lambda = optimal_lambda,
                 initial.estimates = initial.estimates,
                 lambda.initial = lambda.initial,
                 alpha = alpha,
                 beta.initial = beta.initial,
                 standardize = standardize,
                 standardize.response = standardize.response,
                 lower.limits = lower.limits,
                 upper.limits = upper.limits)

  # Return regression estimates (no intecept) and optimal lambda
  out <- list(
    beta = fit$beta,
    shrinkageFactors = fit$shrinkageFactors,
    lambda = optimal_lambda
  )
  return(out)
}

#' @title Estimates the value of lambda for NNG based on AIC and BIC
#'
#' @description
#' The value of lambda that yields the smallest AIC or BIC is considered the best lambda
#' for tuning the nonnegative garrote.
#' @inheritParams nngfit
#' @param lambda.initial The tuning parameter for the initial estimates using
#' lasso or ridge regularization. If \code{initial.estimates = "default"},
#' `lambda.initial` can be NULL since OLS estimates do not require a tuning
#' parameter. However, if \code{initial.estimates = "ridge"} or
#' \code{initial.estimates = "lasso"}, the tuning parameters must be supplied.
#' @param beta.initial A vector of initial estimates obtained from standardized
#' predictor matrix x. The default value is \code{beta.initial = NULL}, and the
#' program will estimate them based on the selected `initial.estimate` and
#' `lambda.initial`.
#'@export
nng_aic_bic<- function(x, y, beta.initial, lambda.initial, criterion,
                       initial.estimates, lower.limits, upper.limits, lambda,
                       alpha, sigma2){

  # set up data
  x <- as.matrix(x)
  y <- as.numeric(y)

  # sample size
  n <- as.integer(dim(x)[1L])

  # product of initial estimates and x. use concept of sparse matrices for speed
  matb <- Matrix::Matrix(diag(beta.initial,nrow = length(beta.initial)), sparse = TRUE)

  x.tilde <- Matrix::crossprod(Matrix::Matrix(Matrix::t(x), sparse = TRUE), matb)

  # Obtain the shrinkage factors for each lambda using glmnet() function.
  fit <- glmnet::glmnet(x = x.tilde,
                        y = y,
                        standardize = FALSE,
                        intercept = FALSE,
                        lambda = lambda,
                        lower.limits = lower.limits, # Important: c>=0
                        family = "gaussian",
                        weights = NULL,
                        offset = NULL,
                        alpha = 1
                        )

  # Matrix of shrinkage factors and corresponding
  cnew <- fit$beta

  # linear predictor
  xb <- as.matrix(Matrix::crossprod(Matrix::t(x.tilde), cnew))

  # residual for each value of lambda
  #res = as.matrix(-1*sweep(x = xb, MARGIN = 1, STATS = y, FUN = "-"))
  res <- matrix(y, nrow = nrow(xb), ncol = ncol(xb), byrow = F) - xb

  # residual sum of squares
  rss <- colSums(res ^ 2, na.rm = T)

  # degrees of freedom, no intercept since x and y were centered.
  #df <- fit$df
  df <- calculate_df(x = x, betas = as.matrix(cnew))

  # AIC and BIC calculation
  pn <- (rss / (n * sigma2))
  aic <- pn + (2 / n) * df
  bic <- pn + (log(n) / n) * df

  # Find lambda that yields smallest AIC or BIC
  lambda.aic <- lambda[which.min(aic)]
  lambda.bic <- lambda[which.min(bic)]

  # fit garrote model
  fitm <- switch(criterion,
               "aic" = garrote(x = x, y = y,
                               lambda = lambda.aic,
                               initial.estimates = initial.estimates,# does not play any role since beta.initial is available
                               lambda.initial = lambda.initial,
                               beta.initial = beta.initial,
                               alpha = alpha,
                               standardize = FALSE,
                               lower.limits = lower.limits,
                               upper.limits = upper.limits
                               ),
               "bic" = garrote(x = x, y = y,
                               lambda = lambda.bic,
                               initial.estimates = initial.estimates,
                               lambda.initial = lambda.initial,
                               beta.initial = beta.initial,
                               alpha = alpha,
                               standardize=FALSE,
                               lower.limits = lower.limits,
                               upper.limits = upper.limits
                               )
               )

  # return parameters of interest
  parms <-list(
      beta = fitm$beta,
      shrinkageFactors = fitm$shrinkageFactors,
      lambda = fitm$lambda,
      aic = min(aic),
      bic = min(bic)
    )
  return(parms)
}


#' @title Calculate the values of lambdas for the nonnegative garrote
#'
#' @description
#' This function calculates the sequence of lambdas with a specified length, nlambda. It
#' assumes that both x and y are already standardized.
#' @param x A standardized matrix of predictors.
#' @param y A standardized vector of the response variable.
#' @param beta.initial Initial estimates for the standardized matrix x.
#' @param nlambda The number of lambdas to be returned.
#' @param epsilon The smallest value used in calculating the minimum lambda. The default value is
#' 0.0001, referred to elsewhere as `lambda.min.ratio`.
#' @param alpha The mixing parameter, where alpha = 1 represents the usual nonnegative
#' garrote, while alpha = 0 represents the garrote with a ridge-like penalty.
#' @return Returns a sequence of lambdas with a length of nlambda.
#' @importFrom Matrix crossprod
#' @export
lambda_seq_garrote <- function(x, y, beta.initial,nlambda = 100,
                              epsilon =  0.0001, alpha = 1){


  # Calculate new matrix of x
  matb <- Matrix::Matrix(diag(beta.initial,nrow = length(beta.initial)), sparse = TRUE)

  xnew <- Matrix::crossprod(Matrix(Matrix::t(x), sparse = TRUE), matb)

  # calculate lambda-max that gives all estimates equal to 0
  #lambda_max = max(abs(t(xnew)%*%y))/length(y)
  #lambda_max <- max(abs(colSums(as.matrix(xnew) * y, na.rm = T))) / length(y)
  alpha <- ifelse(alpha == 0, 0.001, alpha)
  lambda_max <- max(abs(Matrix::crossprod(xnew, y)), na.rm = TRUE) / (length(y) * alpha)

  # calculate lambda-min that gives all estimates equal to nonzero
  lambda_min <- epsilon * lambda_max

  # choose the grid of lambdas to be equidistant on the log-scale
  lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

  # return lambdas
  return(lambda)
}

#' Extract coefficients from nngfit object
#'
#' @param object Fitted \code{"nngfit"} object
#'@param ... not used at the moment
#' @method coef nngfit
#' @export
coef.nngfit <- function(object, ...) {
  object$beta
}

#' make predictions from "nngfit" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"nngfit"} object.
#'
#' @param object Fitted \code{"nngfit"} model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict nngfit
#' @export
predict.nngfit <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  newx %*% coef.nngfit(object)
}

#' @title Calculates a sequence of lambdas for adaptive lasso regularization.
#'
#' @description Calculates the sequence of lambdas of length `nlambda`. If gamma
#' is a vector, it calculates all the lambdas for each gamma, then finds the
#' overall maximum and minimum, and uses them to compute the final sequence
#' of lambdas on a logarithmic scale. It assumes that x and y are standardized.
#' @param x A standardized matrix of predictors
#' @param y A standardized vector of response
#' @param nlambda Number of lambdas
#' @param epsilon The ratio of smallest value of the generated lambda
#' sequence (say lambda.min) to lambda.max. See lambda.min.ratio
#' in \code{[cv.glmnet]} for more details
#' @param beta.initial Initial estimates for the standardized matrix x.
#' @param gamma The tuning parameter for adaptive lasso weights. The default
#' value is NULL, and the program uses seq(0.5, 2, 0.5).
#' @param ignore.gamma Specifies whether to ignore values of gamma supplied by the
#' user and use default values seq(0.5, 2, 0.5)
#' @return Returns sequence of lambdas
#' @export
lambda_seq_alasso <- function(x, y, beta.initial,nlambda = 100,
  epsilon =  0.0001, gamma = NULL, ignore.gamma = FALSE){

  # Set default gammas
  if (is.null(gamma)) {
    gamma <- seq.int(0.5, 2, 0.5)
  }

  # Ignore the sequence of gamma provided
  if (ignore.gamma) {
    gamma <- seq.int(0.5, 2, 0.5)
  }


  # number of gammas
  k <- length(gamma)
  ny <- length(y)

  # transpose x and make it sparse
  x <- Matrix::Matrix(Matrix::t(x),sparse = TRUE)

  lambda_list <- vector(mode = "list", length = k)

  for (i in 1:k) {

    # if all initial estimates are 0, error will be return by seq()
    ## Calculate new matrix of x
    matb <- Matrix::Matrix(diag(abs(beta.initial)^gamma[i]), sparse = TRUE)

    xnew <- Matrix::crossprod(x,matb)

    # calculate lambda-max that gives all estimates equal to 0
    #lambda_max <- max(abs(t(as.matrix(xnew))%*%y))/length(y)-little bit slower
    #lambda_max <- max(abs(colSums(as.matrix(xnew) * y)), na.rm = T) / length(y)
    lambda_max <- max(abs(Matrix::crossprod(xnew, y)), na.rm = TRUE) / ny

    # calculate lambda-min that gives all estimates equal to nonzero
    lambda_min <- epsilon * lambda_max

    # choose the grid of lambdas to be equidistant on the log-scale
    lambda_list[[i]] <- exp(seq(log(lambda_max), log(lambda_min),
      length.out = nlambda))
  }

  # Overall maximum lambda
  lambdamax <- max(unlist(lambda_list, use.names = FALSE))

  # Overall minimum
  lambdamin <- min(unlist(lambda_list, use.names = FALSE))

  # choose the grid of lambdas to be equidistant on the log-scale
  # this will return error when all initial estimates are 0
  lambda <- exp(seq(log(lambdamax), log(lambdamin), length.out = nlambda))

  # return lambdas
  return(lambda)
}

