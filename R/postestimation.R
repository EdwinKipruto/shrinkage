#' @title Estimates parameter-wise shrinkage factors using k-fold
#' cross-validation. 
#' @description
#' This function is a wrapper for pwfs(), which is used to estimate 
#' parameter-wise shrinkage factors based on the approach proposed by 
#' Sauerbrei (1999) and Kipruto and Sauerbei (2024).
#' @param x A standardized matrix of predictors with dimensions nobs x nvars, where 
#' each row represents an observation vector. If not standardized, set standardize to
#' TRUE.
#' @param y A numeric, centered quantitative vector of responses. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. 
#' @param standardize Specifies whether to standardize the predictors matrix 
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. This standardization step can be useful to ensure 
#'  that the predictors are on a comparable scale. By default, 
#'  \code{standardize=FALSE}, indicating that no standardization will be performed.
#'  This assumes that users have already standardized their data.
#' @param varx The names of variables selected by a selection procedure. The 
#' default value is NULL, and all column names of x are used.
#' @param  nfolds The number of cross-validation. if nfolds = n it returns
#'  LOOCV. Default is 10-fold CV
#' @param foldid An optional vector of values between 1 and nfolds identifying
#'  what fold each observation is in while conducting cross-validation. 
#' @param nonnegative Specify whether nonnegative shrinkage factors should be 
#' estimated. Default is FALSE, meaning that negative shrinkage factors
#' might be estimated.
#' @return Returns a list with the following components:
#' \item{ShrunkenRegCoef:}{Shrunken regression estimates for the  variables `varx` without an
#' intercept.}
#' \item{shrinkagefactors:}{Shrinkage factors for the variables `varx`.}
#' @export
shrinkcoef <- function(x, y, 
                       varx = NULL,
                       nfolds = 10, 
                       foldid = NULL,
                       standardize = FALSE, 
                       nonnegative = FALSE) {
  
  # set up data
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # x must have column names
  xnames <- colnames(x)
  
  if (is.null(xnames)) {
    stop("x must have column names.")
  }
  
  # if varx is null use the column names of x
  if (is.null(varx)) {
    varx <- xnames
  }
  
  # subset x using the names of the selected variables denoted by varx
  x <- x[, varx, drop = FALSE]

  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and means of x
    sdx <- apply(x, 2, sd_x)
    
    colx_means <- colMeans(x, na.rm = TRUE)
    
    # center x by means and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)
    
    # center the response variable y
    y <- y - mean.default(y, na.rm = TRUE)
  }
  
  # sample size (n) and number of variables (nvars)
  np <- dim(x)
  n <- as.integer(np[1])
  nvars <- as.integer(np[2])
  
  # Split the data into approximately equal groups
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }
  
  # Important for LOOCV when foldid is provided
  nfolds <-  length(unique(foldid))

  # save partial linear predictors and y since data are sampled randomly
  xb <- ycv <- vector(mode = "list", length = nfolds)

  for (i in 1:nfolds) {
    # Fit a linear model to the nfolds-1 training samples
    coefx <- .lm.fit(x = x[foldid != i, , drop = FALSE],
                     y = y[foldid != i])$coefficients
    
    # obtain partial linear predictors of the ith fold
    #xb[[i]] <- x[foldid==i,,drop=F]%*%diag(coefx, nrow = dim(x)[2])
    xb[[i]] <- crossprod(t(x[foldid == i, , drop = FALSE]), 
                         diag(coefx, nrow = length(coefx))
                         )
    # save y of the ith fold
    ycv[[i]] <- y[foldid == i]
  }
  
  # combine the partial predictors into a single matrix of dimension n * nvars
  predmatrix <- as.matrix(do.call(rbind, xb))
  colnames(predmatrix) <- varx
  
  y2 <- unlist(ycv, use.names = FALSE)
  
  # Estimate nonnegative shrinkage factors: We prefer glmnet due to its speed,
  # but it does not support one variable, which can occur in subset selection.
  # Therefore, for one variable, we will use nnls() from the nnls package. The
  # results from nnls() are almost identical to those from glmnet.
  if (nonnegative) {
    if (ncol(predmatrix) == 1) {
      # Lawson-Hanson NNLS implemention of non-negative least squares
      shrinkagefactors <- setNames(as.numeric(nnls::nnls(A = predmatrix,
                                                         b = y2)$x),
                                   varx)
      
    } else {
      # we can estimate shrinkage factors direct using glmnet
      shrinkagefactors <- glmnet::glmnet(x = predmatrix,
                                         y = y2,
                                         lambda = 0,
                                         lower.limits = 0,
                                         family = "gaussian",
                                         standardize = F)$beta
      
      shrinkagefactors <- setNames(as.numeric(shrinkagefactors), varx)
    }
  } else{
    # original PWSF. This approach can result in undesirable -ve shrinkage factors.
    shrinkagefactors <-setNames(.lm.fit(x = predmatrix,
                                        y = y2)$coefficients,
                                varx)
    
  }
  # unshrunken regression estimates for the selected variables
  coef_selected_model <- .lm.fit(x = x, y = y)$coefficients
  
  # shrunken regression estimates in standardized scale without intercept
  ShrunkenRegCoef <- setNames(coef_selected_model * shrinkagefactors, varx)
  
  # return
  fit <- list(ShrunkenRegCoef = ShrunkenRegCoef, # No intercept
              shrinkageFactors = shrinkagefactors
              )
  # assign class for use in generic functions
  class(fit) = "shrinkcoef"
  
  return(fit)
}

#' @title Estimate parameter-wise shrinkage factors for variables
#' @description Estimates parameter-wise shrinkage factors and
#' corresponding shrunken regression estimates using k-fold cross-validation (CV)
#' or leave-one-out CV. Parameter-wise shrinkage (PWS) approach is an extension
#' of the global shrinkage proposed by Sauerbrei (1999). Unlike global
#' shrinkage, it shrinks regression coefficients differently and was originally
#' proposed for use after model selection. However, Kipruto and Sauerbrei
#' (2024) modified the approach by imposing non-negativity constraints on the 
#' shrinkage factors to allow its usage in both full and selected models.
#' @param x A standardized matrix of predictors with dimensions nobs x nvars, where 
#' each row represents an observation vector. If not standardized, set standardize to
#' TRUE.
#' @param y A numeric, centered quantitative vector of responses. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. 
#' @param standardize Specifies whether to standardize the predictors matrix 
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#'  standardized. Additionally, the response variable \code{`y`} will be centered
#'  by subtracting its mean. This standardization step can be useful to ensure 
#'  that the predictors are on a comparable scale. By default, 
#'  \code{standardize=FALSE}, indicating that no standardization will be performed.
#'  This assumes that users have already standardized their data.
#' @param varx Names of the variables selected by a procedure. Used to subset 
#' the predictor matrix \code{`x`}. The default value is NULL, and all column 
#' names of x are used.
#' @param  nfolds The number of cross-validation. if nfolds = n it returns
#'  LOOCV. The default value is \code{nfolds = 10}, which conducts 10-fold CV.
#' @param foldid An optional vector of values between 1 and nfolds identifying
#'  what fold each observation is in while conducting cross-validation. 
#' @param nonnegative Specify whether nonnegative shrinkage factors should be 
#' estimated. Default is FALSE, meaning that negative shrinkage factors might be 
#' estimated.
#' @param loocv Specifies whether leave-one-out cross-validation (LOOCV) should be
#' used to estimate shrinkage factors. The default is set to FALSE, indicating that
#' k-fold cross-validation is performed.
#' @return Returns a list with the following components:
#'\item{ShrunkenRegCoef}{Shrunken regression estimates for the selected variables,
#' with coefficients of zero assigned to the unselected variables.}
#'\item{ShrinkageFactors}{Shrinkage factors for the selected variables, with 
#'shrinkage factors of zero assigned to the unselected variables.}
#'\item{variables_selected}{The variables chosen for estimating shrinkage
#' factors (varx)}
#'\item{nvar}{The number of variables selected by the procedure (length of 
#'\code{"varx"}).}
#'\item{x}{The standardized design matrix of all predictor variables.}
#'\item{y}{The centered response variable.}
#' @export
pwsf <- function(x, y,
                 standardize = FALSE,
                 varx = NULL,
                 nfolds = 10,
                 foldid = NULL,
                 nonnegative = FALSE,
                 loocv = FALSE
                 ){
  
  # Set up data
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # assert that x must have names
  xnames <- colnames(x)
  
  if (is.null(xnames)) {
    stop("x must have column names")
  }
  
  # Work on selected x variables, which can also be variables from the full model.
  if (is.null(varx)) {
    varx <- xnames
  }
  
  # sample size and number of variables
  dimx <- dim(x)
  n <- as.integer(dimx[1L])
  nvar <- as.integer(dimx[2L])
  
  # standardize x and center y if necessary
  if (standardize) {
    # standard deviations and column means of x
    sdx <- apply(x, 2, sd_x)
    
    colx_means <- colMeans(x, na.rm = TRUE)
    
    # center x and scale x by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)
    
    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }
  
  # Deal with foldid for splitting data in cross-validation
  if (is.null(foldid)) {
    if (loocv) {
      nfolds <- n
    }
    foldid <- sample(rep(seq(nfolds), length = n))
  } else {
    # we will overwrite foldid if loocv = T
    if (loocv) {
      nfolds <- n
      foldid <- sample(rep(seq(nfolds), length = n))
    } else {
      foldid <- foldid
    }
  }
  
  # Estimate shrinkage factors and shrunken regression estimates-no intercept
  fitx <-  shrinkcoef(x = x,
                      y = y, 
                      varx = varx, 
                      foldid = foldid,
                      nonnegative = nonnegative, 
                      standardize = FALSE
                      )
  
  # return regression estimates like penalized methods: padded zeros
  index <- match(varx, xnames)
  betas <- setNames(rep(0, nvar), xnames)
  betas[index] <- fitx$ShrunkenRegCoef
  
  # return shrinkage factors like penalized methods: zeros means variable not selected
  shrinkagefactors <- setNames(rep(0, nvar), xnames)
  shrinkagefactors[index] <- fitx$shrinkageFactors

  fit <- list(
    ShrunkenRegCoef = betas, # no intercept
    shrinkageFactors = shrinkagefactors,
    variables_selected = varx,
    nvar = length(varx),
    x = x,
    y = y
  )
  
  # assign class
  class(fit) = "pwsf"
  
  return(fit)
}

#' Extract shrunken regression estimates from pwsf object
#'
#' @param object Fitted \code{"pwsf"} object
#' @param ... not used at the moment
#' @method coef pwsf
#' @export 
coef.pwsf <- function(object, ...) {
  object$ShrunkenRegCoef
}

#' make predictions from a "pwsf" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"pwsf"} object.
#' @param object Fitted \code{"pwsf"} object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict pwsf
#' @export
predict.pwsf <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  newx %*% coef.pwsf(object)
}

#' @title Estimate global shrinkage factors
#' @description Estimates global shrinkage factors and calculates corresponding 
#' shrunken regression estimates using either k-fold cross-validation (CV) or 
#' leave-one-out CV. Global shrinkage involves shrinking all regression coefficients
#' by the same factor, as proposed by van Houwelingen and le Cessie (1990) and 
#' Verweij and van Houwelingen (1993).
#' @param x A standardized matrix of predictors with dimensions nobs x nvars, where 
#' each row represents an observation vector. If not standardized, set standardize to
#' TRUE.
#' @param y A numeric, centered quantitative vector of responses. If not centered,
#' set `standardize` to TRUE; this will standardize x but center y. 
#' @param standardize Specifies whether to standardize the predictors matrix 
#' \code{`x`} to have mean 0 and unit variance. If set to TRUE, \code{`x`} will be
#' standardized. Additionally, the response variable \code{`y`} will be centered
#' by subtracting its mean. This standardization step can be useful to ensure 
#' that the predictors are on a comparable scale. By default, 
#' \code{standardize=FALSE}, indicating that no standardization will be performed.
#' This assumes that users have already standardized their data.
#' @param varx Names of the variables selected by a procedure. Used to subset 
#' the predictor matrix \code{`x`}. The default value is NULL, and all column
#' names of x are used.
#' @param nfolds Number of cross-validation folds. The default is 10; setting 
#' `nfolds = n` corresponds to leave-one-out cross-validation.
#' @param foldid An optional vector of values between 1 and \code{"nfolds"} that 
#' identifies the fold to which each observation belongs. The default is 
#' NULL, and the program will generate it.
#' @param nonnegative Specify whether nonnegative shrinkage factors should be 
#' estimated. Default is FALSE, meaning that negative shrinkage factors might
#' be estimated. Negative shrinkage factor is hardly estimated unless in extreme
#' situations where SNR is extremely low.
#' @param nonnegative Specify whether nonnegative shrinkage factors should be 
#' estimated. The default is FALSE, allowing for the estimation of negative 
#' shrinkage factors. Negative shrinkage factors are rarely estimated unless 
#' in extreme situations where Signal-to-noise ratio is extremely low.
#' @return Returns a list with the following components:
#'\item{ShrunkenRegCoef}{Shrunken regression estimates for the selected variables,
#' with coefficients of zero assigned to the unselected variables.}
#'\item{ShrinkageFactors}{Shrinkage factor for the selected variables, with 
#'shrinkage factors of zero assigned to the unselected variables.}
#'\item{variables_selected}{The variables chosen for estimating shrinkage
#' factors (varx)}
#'\item{nvar}{The number of variables selected by the procedure (length of 
#'\code{"varx"}).}
#' @export
global <- function(x, y,
                   varx = NULL,
                   nfolds = 10,
                   foldid = NULL,
                   standardize = FALSE,
                   nonnegative = FALSE){
  
  #set up data
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # number of all variables in x before subsetting
  nv <- dim(x)[2L]
  
  # assert that x must have names before subsetting
  xnames <- colnames(x)
  
  if (is.null(xnames)) {
    stop("x must have column names")
  }
  
  # work only on the selected x variables. varx can also be from a full model
  if (is.null(varx)) {
    varx <- xnames
  }
  
  # subset x based on varx
  x <- x[, varx, drop = FALSE]
  
  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of matrix x
    sdx <- apply(x, 2, sd_x)
    
    colx_means <- colMeans(x, na.rm = T)
    
    # center x by mean and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)
    
    # center the response y by mean
    y <- y - mean.default(y, na.rm = TRUE)
  }
  
  # sample size and number of variables after subsetting using varx
  np <- dim(x)
  n <- as.integer(np[1])
  nvars <- as.integer(np[2])
  
  # Split the data into approximately equal groups for cv
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = n))
  }
  
  # Important for LOOCV when foldid is provided
  nfolds <- length(unique(foldid))
  
  # save partial linear predictors and y since data are sampled randomly
  xb <- ycv <- vector(mode= "list", length = nfolds)

  for (i in 1:nfolds) {
    # Fit a linear model to the k-1 training samples
    coefx <- .lm.fit(x = x[foldid != i, , drop = FALSE],
                     y = y[foldid != i])$coefficients
    
    # find the product of the estimates with x of the ith fold
    xb[[i]] <- crossprod(t(x[foldid == i, , drop = FALSE]), coefx)
    ycv[[i]] <- y[foldid == i]
  }
  
  # combine the linear predictors into a single matrix of dim = n * 1
  predmatrix <- as.matrix(do.call(rbind, xb))
  colnames(predmatrix) <- "eta"
  y2 <- unlist(ycv, use.names = FALSE)
  
  # regress y2 against predmatrix to obtain global shrinkage factor
  if (nonnegative) {
    shrinkagefactor <- as.numeric(nnls::nnls(A = predmatrix, 
                                              b = y2)$x)
  } else {
    shrinkagefactor <- as.numeric(.lm.fit(x = predmatrix, 
                                           y = y2)$coefficients)
  }

  # unshrunken standardized regression estimated for the selected variables
  regcoef <- .lm.fit(x = x, y = y)$coefficients
  
  # shrunken regression estimates in standardized scale
  ShrunkenRegCoef <-  regcoef*shrinkagefactor
  
  # Return the betas and shrinkage factors in a similar manner to penalized 
  # methods: the betas of eliminated variables are padded with zeros.
  index <- match(varx, xnames)
  
  # shrunken regression estimates for all variables in x
  ShrunkenRegCoef_all_var <- setNames(rep(0, nv), xnames)
  ShrunkenRegCoef_all_var[index] <- ShrunkenRegCoef
  
  # shrinkage factors for all variables in x
  shrinkagefactors_all_var <- setNames(rep(0, nv), xnames)
  shrinkagefactors_all_var[index] <- shrinkagefactor
  
  # return
  fit <- list(ShrunkenRegCoef = ShrunkenRegCoef_all_var, # no intercept
              shrinkageFactors = shrinkagefactors_all_var,
             variables_selected = varx,
             nvar = length(varx)
             )
  
  # assign class
  class(fit) = "global"
  
  return(fit)
}


#' @title Estimates parameter-wise shrinkage factors using Breiman (1995) method
#' 
#' @description This method uses a quadratic penalty term on the shrinkage factors
#' (Breiman 1995). The original approach did not impose non-negativity constraints on 
#' the shrinkage factors. In this implementation, both options are available through
#' the "lower.limits" parameter.
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
#' standardized to have unit variance. This option divides \code{`y`} by its standard
#' deviation.
#' @param lambda  Optional vector of tuning parameter. Default is NULL and 
#' the program will generate its sequence
#' @param foldid A vector of values between 1 and nfolds identifying what fold
#' each observation is in while conducting cross-validation. Default is NULL
#' and the program will generate it.
#' @param nfolds Number of folds - default is 10. If k = n then k-fold is equal to
#' leave-one-out cv.
#' @param  varx Variables selected by a selection procedure. The default value 
#' is NULL, and all column names of \code{x} are used.
#' @param lambda.initial The tuning parameter for the initial estimates using
#' lasso or ridge regularization. If a vector of lambdas are provided, the
#' optimal lambda for initial estimates will be returned.
#' @param beta.initial A vector of initial estimates obtained from standardized 
#' predictor matrix x. The default value is \code{beta.initial = NULL}, and the 
#' program will estimate them based on the chosen `initial.estimate` and 
#' `lambda.initial`. If supplied, the number of nonzero elements must be at least 2.
#' @param initial.estimates Specifies the type of initial estimates to use for
#' constructing weights. The \code{initial.estimates = "default"}
#' option uses \code{"OLS"} estimates. Other options include 
#' \code{initial.estimates ="ridge"} or \code{initial.estimates ="lasso"} 
#' estimates, which are calculated using the user-supplied \code{`lambda.initial`.}
#' @param lower.limits Vector of lower limits for each shrinkage factor. 
#' The default is 0, and nonnegative shrinkage factors will be estimated.
#' @param upper.limits Vector of upper limits for the shrinkage factors. 
#' The default is Inf, and unbounded shrinkage factors will be estimated. 
#' @param type.measure Loss function to use for cross-validation. Currently, 
#' two options are available. The default option is \code{type.measure = "mse"},
#' which corresponds to squared error. Another option is \code{type.measure = "mae"}
#' (mean absolute error).
#' @param nlambda The number of lambda values; default is 100
#' @param lambda.min.ratio lambda can be provided if the user wants to specify 
#' the lambda sequence, but typical usage is for the program to construct the 
#' lambda sequence on its own. When automatically generated, the lambda sequence
#' is determined by lambda.max and lambda.min.ratio. The latter is the ratio of
#' smallest value of the generated lambda sequence (say lambda.min) to lambda.max.
#' The program generates nlambda values linear on the log scale from lambda.max 
#' down to lambda.min. lambda.max is not user-specified but is computed from the
#' input x and y: it is the smallest value for lambda such that all the 
#' coefficients are zero. Default is lambda.min.ratio = 0.00000001
#' @return Returns a list with the following components:
#' \item{ShrunkenRegCoef}{Estimated shrunken coefficients}
#' \item{nvar}{The number of variables selected}
#' \item{shrinkagefactors}{ Shrinkage factors for the selected variables}
#' \item{x}{Design matrix of the selected variables}
#' @export
breiman <- function(x, y, 
                    varx = NULL,
                    lambda = NULL, 
                    nfolds = 10, 
                    foldid = NULL, 
                    nlambda = 100,
                    beta.initial = NULL,
                    lambda.initial= NULL,
                    initial.estimates = c("default", "ridge", "lasso"),
                    lambda.min.ratio = 1e-08,
                    type.measure = c("mse", "mae"),
                    lower.limits = 0,
                    upper.limits = Inf, 
                    standardize = FALSE,
                    standardize.response = FALSE){
  
  # match arguments
  initial.estimates <- match.arg(initial.estimates)
  type.measure <- match.arg(type.measure)
  
  # beta.initial must have at least two nonzero coefficients
  if (!is.null(beta.initial) && sum(beta.initial!=0)<2)
    stop("The number of nonzero elements in `beta.initial` must be at least 2")
  
  # set up data
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  # number of all variables in x before subsetting
  nv <- dim(x)[2L]
  
  # assert that x must have names before subsetting
  xnames <- colnames(x)
  
  if (is.null(xnames))
    stop("x must have column names")
  
  # work only on the selected x variables. varx can also be from a full model
  if (is.null(varx)) {
    varx <- xnames
  }
  
  # subset x using varx
  x <- as.matrix(x[, varx, drop = FALSE])
  
  # standardize x and center y if necessary
  if (standardize) {
    # standard deviation and mean of matrix x
    sdx <- apply(x, 2, sd_x)
    colx_means <- colMeans(x, na.rm = TRUE)
    
    # center x by mean and scale by standard deviation
    x <- sweep(x, 2L, colx_means, "-", check.margin = FALSE)
    x <- sweep(x, 2L, sdx, `/`, check.margin = FALSE)
    
    # center the response y
    y <- y - mean.default(y, na.rm = TRUE)
  }
  
  # Standardize the response variable to have sd of 1
  if(standardize.response){
    y <- y / sd_x(y)
  }
  
  # estimate lambda via cross-validation
  cvout <- cv.garrote(x = x, y = y, 
                      nfolds = nfolds,
                      foldid = foldid, 
                      lambda = lambda,
                      nlambda = nlambda,
                      lambda.initial = lambda.initial,
                      beta.initial = beta.initial,
                      alpha = 0,
                      type.measure = type.measure,
                      initial.estimates = initial.estimates,
                      lower.limits = lower.limits,
                      upper.limits = upper.limits,
                      standardize = FALSE,
                      standardize.response = FALSE,
                      lambda.min.ratio =lambda.min.ratio)
  
  # use optimal lambda to estimate the shrinkage factors
  fit <- garrote(x = x, y = y,
                 alpha = 0,
                 lambda = cvout$lambda.min,
                 lambda.initial = cvout$lambda.initial,
                 initial.estimates = initial.estimates,
                 lower.limits = lower.limits,
                 upper.limits = upper.limits,
                 standardize = FALSE,
                 standardize.response = FALSE)
  
  # return betas like penalized methods: padded zeros

  # match the selected variables and all variables
  index <- match(varx, xnames)
  
  # shrunken regression estimates
  ShrunkenRegCoef <- setNames(rep(0, nv), xnames)
  ShrunkenRegCoef[index] <- fit$beta
  
  # shrinkage factors
  shrinkagefactors <- setNames(rep(0, nv), xnames)
  shrinkagefactors[index] <- fit$shrinkageFactors
  
  
  fit <- list(ShrunkenRegCoef = ShrunkenRegCoef, # regression estimates without intercept
              shrinkageFactors = shrinkagefactors,
              nvar = length(varx),
              x = x, 
              y = y
              )
  
  # assign class to use generic functions
  class(fit) = "breiman"
  return(fit)
}

#' Extract coefficients from a breiman object
#'
#' @method coef breiman
#' @param object Fitted \code{"breiman"} object
#' @param ... not useful at the moment
#' @export
coef.breiman <- function(object, ...) {
  object$ShrunkenRegCoef
  
}

#' make predictions from a "breiman" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' from a fitted \code{"breiman"} object.
#' @param object Fitted \code{"breiman"} object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix.
#' @param ... not used at the moment
#' @method predict breiman
#' @export
predict.breiman <- function(object, newx = NULL, ...) {
  if (is.null(newx)) {
    newx <- object$x
  }
  newx %*% coef.breiman(object)
}

#' Calculate shrinkage factors between two beta matrices
#'
#' This function takes two matrices of beta coefficients and calculates the 
#' shrinkage factors by dividing beta1 by beta2. It can handle cases where 
#' certain variables need to be excluded in the calculation.
#'
#' @param beta1 The beta matrix used to calculate shrinkage factors (e.g., full
#' model or backward elimination without post-estimation shrinkage).
#' @param beta2 The beta matrix from a regression procedure that has already 
#' employed shrinkage (e.g., backward elimination with global shrinkage).
#' @param exclude A numeric vector specifying the column indices of variables to
#' be excluded in both beta matrices.
#'
#' @return A matrix of shrinkage factors with the same dimensions as the input
#' matrices, or a subset of it if variables are excluded.
#'
#' @examples
#' beta1 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
#' beta2 <- matrix(c(0.5, 0.5, 2, 2, 3, 3), nrow = 2)
#' shrink_betas(beta1, beta2)
#' @export
shrink_betas <- function(beta1, beta2, exclude = NULL) {
  
  # set up data
  beta1 <- as.matrix(beta1)
  beta2 <- as.matrix(beta2)
  
  if (!is.null(exclude)) {
    if (!is.numeric(exclude)) {
      stop("exclude must be a numeric vector.")
    }
    if (any(exclude < 1) || any(exclude > ncol(beta1))) {
      stop("Invalid exclude column indices.")
    }
  }
  
  outx <- beta1 / beta2
  outx[is.nan(outx)] <- 0
  
  if (!is.null(exclude)) {
    outx <- outx[, -exclude, drop = FALSE]
  }
  
  return(outx)
}


