#' @title Generates simulation data.
#' @description Generate a predictor matrix \code{x} and response vector \code{y} based 
#' on a specified setup. Two sets of predictor-response pairs are generated: one for 
#' training and one for testing. The training \code{x} is standardized to have a mean of 0 
#' and unit variance, while the training response is centered and standardized if necessary. 
#' The test data are then standardized using statistics derived from the training data.
#' @param n The number of training observations in the training data (sample size). Default is 100
#' @param ntest The number of observations in the testing data. Default is 100,000 observation
#' @param nvar The number of predictors in training and testing dataset. Default is 15
#' @param corrtype The type of correlation matrix used. Four options are available:\cr
#' \code{c1}: Correlation matrix from Houwelingen and Sauerbrei (2013)\cr
#' \code{c2}: AR1(0.3) autocorrelation\cr
#' \code{c3}: AR1(0.8) autocorrelation\cr
#' \code{c4}: Correlation from educational Body fat dataset
#' \code{c5}: Block correlation matrix\cr
#' Note that for AR1 autocorrelation, the correlation between the ith and jth
#' variables is defined as AR1(rho) = rho^|i-j|.
#' @param withinSignalCorr,withinNoiseCorr,betweenSignalNoiseCorr parameters of block
#' correlation type (c5). See corr_type() for details. description
#' @param betatype Specifies the type of true regression coefficients. Please refer
#' to the simulation protocol for more detailed information. The nonzero coefficients
#' in the regression model are fixed at a constant value of 7.
#' @param standardize.response Specifies whether the training response 
#' variable should be standardized. Default is FALSE and y is only centered by mean
#' @param snr Specifies the desired signal-to-noise ratio (SNR), i.e., 
#' \eqn{SNR = \frac{var(mu)}{\sigma^2}} where \eqn{var(mu)} is the variance
#' of linear predictor and \eqn{\sigma^2} is the error variance. The error
#' variance is set so that the given SNR is achieved. Default is SNR = 1.
#' @return A list with the following components:
#'   \item{\code{x}}{Standardized training predictors}
#'   \item{\code{y}}{Centered (standardized) training response variable} 
#'   \item{\code{xtest}}{Standardized testing predictors}
#'   \item{\code{ytest}}{Centered testing response variable}
#'   \item{\code{Sigma}}{Covariance matrix used to generate the data}
#'   \item{\code{beta}}{True regression coefficients used to generate the data}
#'   \item{\code{sd_y}}{Standard deviation of the training outcome}
#' @details The data model is: \eqn{Y \sim N(X\beta, \sigma^2 I)}. The predictor 
#' variables have covariance matrix `Sigma`
#' @references
#' The simulation setup is based on the study titled "Comparison of variable selection 
#' procedures and investigation of the role of shrinkage: a simulation protocol" by 
#' [Kipruto and Sauerbrei(2022)].
#' @export 
sim_data <- function(n = 100, 
                     nvar = 15, 
                     ntest = 100000, 
                     corrtype = c("c1", "c2", "c3", "c4", "c5"), 
                     betatype = c("a", "b", "c", "d"),
                     withinSignalCorr = 0.8, 
                     withinNoiseCorr = 0.8,
                     betweenSignalNoiseCorr = 0,
                     snr = 1,
                     standardize.response = FALSE) {
  
  # match arguments
  corrtype <- match.arg(corrtype)
  betatype <- match.arg(betatype)
  
  if (nvar < 15) {
    stop("The number of predictors should not be less than 15.")
  }
  # Generate training and validation predictors
  x <- matrix(rnorm(n * nvar), n, nvar)
  xtest <- matrix(rnorm(ntest * nvar), ntest, nvar)
  
  # Introduce correlation to the predictors
  Sigma <- corr_type(nvar = nvar, 
                     type = corrtype, 
                     withinSignalCorr = withinSignalCorr,
                     withinNoiseCorr = withinNoiseCorr,
                     betweenSignalNoiseCorr = betweenSignalNoiseCorr,
                     betatype = betatype)
  obj <- svd(Sigma)
  Sigma.half <- obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
  x <- x %*% Sigma.half
  xtest <- xtest %*% Sigma.half
  
  # mean and standard deviation of training x
  sd_xtrain <- apply(x, 2, sd_x)
  mean_xtrain <- colMeans(x, na.rm = TRUE)
  
  # standardize training x to have mean 0 and unit variance
  x <- sweep(x, 2L, mean_xtrain, "-", check.margin = FALSE)
  x <- sweep(x, 2L, sd_xtrain, `/`, check.margin = FALSE)
  
  # standardize testing x using training x mean and sd
  xtest <- sweep(xtest, 2L, mean_xtrain, "-", check.margin = FALSE)
  xtest <- sweep(xtest, 2L, sd_xtrain, `/`, check.margin = FALSE)
  
  # Assign names to the column of matrices
  colnames(x) = colnames(xtest) = colnames(Sigma) = paste0("x",1:nvar)
  
  # Generate coefficients for the standardized variable x and xtest
  beta <- beta_type(nvar = nvar, type = betatype)
  
  # Set SNR based on sample variance on infinitely large test set
  vmu <- as.numeric(t(beta) %*% Sigma %*% beta)
  sigma <- sqrt(vmu / snr)
  
  # Generate training response variable
  y <- as.numeric(x %*% beta + rnorm(n) * sigma)
  
  # mean and standard deviation of training y
  mean_ytrain <- mean.default(y, na.rm = TRUE)
  sd_y <- sd_x(y)
  
  # Generate testing response variable
  ytest <- as.numeric(xtest %*% beta + rnorm(ntest) * sigma)
  
  # center y train and y test
  y <- y - mean_ytrain
  ytest <- ytest - mean_ytrain
  
  # standardized y train
  if (standardize.response) {
    y <- y / sd_y
  }

  # return data
  enlistx(x, y, xtest, ytest, Sigma, beta, sigma, sd_y)
}

#-------------------------------------------------------------------------------
# Master Function for Conducting Simulation Study
#-------------------------------------------------------------------------------

#' @title Master function for running simulations.
#' @description Run a set of simulations with the specified configuration.
#' @param n The number of observations in the training data(sample size). Default is 100 
#' observations
#' @param ntest The number of observations in the testing data. Default is 100,000 observations
#' @param nvar The number of predictors in the training and testing data.
#' @param corrtype The type of correlation matrix used. Four options are available:
#' \code{c1}: Correlation matrix from Houwelingen and Sauerbrei (2013)\cr
#' \code{c2}: AR1(0.3) autocorrelation\cr
#' \code{c3}: AR1(0.8) autocorrelation\cr
#' \code{c4}: Correlation from educational Body fat dataset
#' \code{c5}: Block correlation structure\cr
#' Note that for AR1 autocorrelation, the correlation between the ith and jth
#' variables is defined as AR1(rho) = rho^|i-j|.
#' @param betatype Specifies the type of true regression coefficients. Please refer
#' to the simulation protocol for more detailed information. The nonzero coefficients
#' in the regression model are fixed at a constant value of 7.
#' @param standardize.response Specifies whether the training response 
#' @param withinSignalCorr,withinNoiseCorr,betweenSignalNoiseCorr parameters of block
#' correlation type (c5). See corr_type() for details.
#' variable should be standardized. Default is FALSE and y is only centered by mean
#' @param reg.funs A list of functions representing the regression
#' procedures to be evaluated in the simulation. Each element of the
#' list must be a function that takes `x`, `y`, `foldid`, and `sigma2` 
#' (i.e., the training predictor matrix, training response vector, vector 
#' of values between 1 and the number of folds identifying what fold each 
#' observation is in while conducting cross-validation, and residual variance 
#' from the full model) as its only four (mandatory) arguments. It must return 
#' an object with associated \code{coef} and \code{predict} methods. The 
#' \code{coef} method must take `object` (the returned object) and return a 
#' vector of coefficients without an intercept. The \code{predict} method must 
#' take `object` and `newx` (the returned object and a new predictor matrix) 
#' and return a matrix of predictions.
#' @param nrep The number of repetitions of which to average the results. 
#' Default is 2000.
#' @param seed The overall random number generator set before repetitions begin 
#' (for the reproducibility of simulation results). The default is 472095.
#' @param verbose Specifies whether intermediate progress should be printed. Default is FALSE.
#' @param file Name of the file to save the simulation results using saveRDS. If
#' set to NULL, no simulation results will be saved. The default value is NULL.
#' @param file.rep The number of repetitions after which intermediate results 
#' are saved. If file.rep is set to 0, it indicates that simulation results 
#' should be saved only at the end with no intermediate saving. The default is 5.
#' @param foldid An optional vector of values between 1 and `nfolds` identifying 
#' what fold each observation is in. If supplied, nfolds can be missing.
#' @param nfolds Number of folds for cross-validation. Default is 10. 
#' Smallest value allowable is nfolds = 3.
#' @param standardize.response Specifies whether the training response variable 
#' should be standardized. If set to TRUE, the model fitted will have both `x` 
#' and `y` standardized, but the coefficients will be back-transformed to the 
#' scale of standardized x.
#' @param snr Specifies the desired signal-to-noise ratio (SNR), i.e., 
#' \eqn{SNR = \frac{var(mu)}{\sigma^2}} where \eqn{var(mu)} is the variance
#' of linear predictor and \eqn{\sigma^2} is the error variance. The error
#' variance is set so that the given SNR is achieved. Default is SNR = 1.
#' @return A list with the following components: TODO
#' \item{err.train}{training error.}
#' \item{err.val}{validation error.}
#' \item{err.test}{test error.} 
#' \item{prop}{test proportion of variance explained}
#' \item{risk}{risk}\item{nzs}{number of selected nonzero coefficients.} \item{fpos}{number of false positives} \item{tpos}{number of true positives}
#' \item{fneg}{number of false negatives.}\item{tneg}{number of true negatives} \item{error.rate}{error rate}
#' \item{opt}{relative optimism (difference in test error and training error, divided by training error).}
#' \item{MCC}{Matthews correlation coefficient.}
#' \item{F1}{F1-score.}
#' Each of length N, where N is the number of regression methods under consideration
#' (the length of reg.funs). The ith element of each list is a matrix of
#' dimension nrep x m, where m = 1 is the number of tuning parameters inherent to
#' the ith method.
#' @seealso \code{\link{sim_data}}
#' @references Kipruto, E. and Sauerbrei, W. (2022). Comparison of variable 
#' selection procedures and investigation of the role of shrinkage in linear
#' regression-protocol of a simulation study in low-dimensional data. Plos one, 
#' 17(10), p.e0271240.
#' @export 
sim_master <- function(n = 100, 
                       nvar = 15, 
                       ntest = 100000, 
                       reg.funs, 
                       nrep = 2000,
                       seed = 472095, 
                       file = NULL, 
                       file.rep = 5,
                       corrtype = c("c1", "c2", "c3", "c4", "c5"), 
                       betatype = c("a", "b", "c", "d"),
                       withinSignalCorr = 0.8, 
                       withinNoiseCorr = 0.8,
                       betweenSignalNoiseCorr = 0,
                       snr = 1,
                       nfolds = 10,
                       foldid = NULL,
                       standardize.response = FALSE,
                       verbose = FALSE) {
  # capture the call
  this.call <- match.call()
  
  # set seed for reproducibility
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  # set parameter for cross-validation grouping for all methods
  if (is.null(foldid)) {
       foldid <- sample(rep(seq(nfolds), length = n))
  }
  # match arguments
  corrtype <- match.arg(corrtype)
  betatype <- match.arg(betatype)
  
  # True regression coefficients
  beta.true <- beta_type(nvar = nvar, type = betatype)
  
  # Number of regression methods
  N <- length(reg.funs)
  
  # Names of regression methods used
  reg.names <- names(reg.funs)
  if (is.null(reg.names)) {
    reg.names <- paste("Method", 1:N)
  }
  
  # Metrics of interest
  err.train = err.val = err.test = nzs = fpos = fneg = tneg = F1 = MCC = fnr = fpr =
  error.rate = risk =  opt = runtime = prop = tpos = rel.test.err = rel.risk = betas = shrinkage = lambda = gamma = vector(mode = "list", length = N)
  
  # Assign names
  names(err.train) = names(err.val) = names(err.test) = names(nzs) = names(risk) =
  names(fpos) = names(tpos) = names(fneg) = names(tneg) = names(F1) = names(MCC) = names(fnr) =
  names(fpr) = names(error.rate) = names(opt) = names(runtime) = names(prop) =
  names(betas) = names(shrinkage) = names(rel.test.err) = names(rel.risk) = 
  names(lambda) = names(gamma) = reg.names
  
  
  for (j in 1:N) {
    err.train[[j]] = err.val[[j]] = err.test[[j]] = nzs[[j]] = fpos[[j]] =
      prop[[j]]=tpos[[j]]= fneg[[j]] = F1[[j]] = MCC[[j]] = risk[[j]] = tneg[[j]] = fnr[[j]] = fpr[[j]]=error.rate[[j]] =
      opt[[j]] = runtime[[j]]  = rel.test.err[[j]] = rel.risk[[j]] = lambda[[j]] = gamma[[j]] = matrix(NA,nrow = nrep,ncol = 1)
    
    betas[[j]]= shrinkage[[j]] = matrix(NA, nrow = nrep, ncol = nvar)
  }
  
  filled = rep(FALSE,N)
  err.null = risk.null = sigma = rep(NA, nrep)
  
  #xx = yy = vector(mode = "list", length = nrep)
  # Loop through the repetitions
  for (i in 1:nrep) {
    if (verbose) {
      cat(sprintf("Simulation %i (of %i) ...\n",i,nrep))
      cat("  Generating data ...\n")
    }
    
    # Generate x, y, xtest and ytest
    simdt <- sim_data(n = n,
                      nvar = nvar,
                      ntest = ntest,
                      corrtype = corrtype,
                      withinSignalCorr = withinSignalCorr,
                      withinNoiseCorr = withinNoiseCorr,
                      betweenSignalNoiseCorr = betweenSignalNoiseCorr,
                      snr = snr,
                      betatype = betatype,
                      standardize.response = standardize.response
                      )
    
    # Training data. x is standardized while y is centered (and standardized if 
    # standardize.response = T)
     x <- simdt$x
     y <- simdt$y
     
     # Test data 
     xtest <- simdt$xtest
     ytest <- simdt$ytest
    
    # Calculate the null risk (betahat = 0):denominator of RR
    risk.null[i] <- diag(t(simdt$beta) %*% simdt$Sigma %*% simdt$beta)
    
    # Numerator of RTE: 
    err.null[i] <- risk.null[i] + simdt$sigma^2
    
    # Residual standard deviation
    sigma[i] <- simdt$sigma
    
    # Loop through the regression methods
    for (j in 1:N) {
      if (verbose) {
        cat(sprintf("Applying regression method %i (of %i) ...\n",j,N))
      }
      
      # Residual variance from full model required for AIC and BIC calculation. 
      fitm <- residual_variance(x = x, 
                                y = y
                                )
      sigma2 <- fitm$sigma2
      beta_full <- fitm$coefficients
      
      # Record system time for each selection method
      tryCatch({
        # Apply the regression method in hand
        runtime[[j]][i] <- system.time({
          # Note: cross-validation is conducted within the reg.funs
          reg.obj <- reg.funs[[j]](x = x,
                                   y = y,
                                   foldid = foldid,
                                   sigma2 = sigma2,
                                   betatypes = betatype
                                  )
        })[1]
        
        # vector of regression estimates without an intercept
        betahat <- coef(reg.obj) 
        
        # Back transform betahat to the original scale of centred y  but standardized x
        # important becauses the beta.true was used with unstandardized y and standardized x
        if (standardize.response) {
          # Sd for x which should be approx 1 because it has already been standardized in sim_data
          sd_xtrain <- apply(x, 2, sd_x)
          
          # The final betahat is for standardized x which we need since xtest is standardized
          betahat <-  (betahat * simdt$sd_y) / sd_xtrain
          
          # Also backtransform beta_full used for estimating shrinkage factors
          beta_full <- (beta_full * simdt$sd_y) / sd_xtrain
        }
        
        # predicted values on the training data
        muhat.train <- predict(reg.obj, newx = x)
        
        # predicted values on the testing data
        muhat.test <- predict(reg.obj, newx = xtest)
        
        # Populate empty matrices for our metrics, of appropriate dimension
        if (!filled[j]) {
          err.train[[j]] = err.val[[j]] = err.test[[j]] = prop[[j]] =risk[[j]]=
            nzs[[j]] = fpos[[j]] =tpos[[j]]= fneg[[j]] = tneg[[j]] = F1[[j]] = MCC[[j]] = 
            fnr[[j]]= fpr[[j]]= error.rate[[j]] = opt[[j]] = rel.test.err[[j]] =
            rel.risk[[j]] = lambda[[j]] = gamma[[j]] = matrix(NA, nrep, 1)
          
          # regression estimates with corresponding shrinkage factors
          betas[[j]] = shrinkage[[j]] = matrix(NA, nrep, nvar)
          filled[j] = TRUE
          # N.B. Filling with NAs is important, because the filled flag could
          # be false for two reasons: i) we are at the first iteration, or ii)
          # we've failed in all previous iters to run the regression method
        }

        # Training data Mean squared error
        err.train[[j]][i, ] <- mean((muhat.train - y) ^ 2, na.rm = TRUE)
        
        # Testing data Mean squared error
        err.val[[j]][i, ] <- mean((muhat.test - ytest) ^ 2, na.rm = TRUE)
        
        # Difference between estimated betas and true betas
        delta <- betahat - simdt$beta
        
        # Risk = Model Error. t(delta) %*% simdt$Sigma %*% delta + beta0^2 if
        # intercept but our intercept = 0
        risk[[j]][i, ] <- t(delta) %*% simdt$Sigma %*% delta
        
        # Test error, compare with error.val
        err.test[[j]][i, ] <- risk[[j]][i, ] + simdt$sigma ^ 2
        
        # proportion of variance explained
        prop[[j]][i,] <- 1-err.test[[j]][i,]/err.null[i]
        
        # Number of nonzero components
        nzs[[j]][i, ] <- sum(betahat != 0)
        
        # number of true positives
        tpos[[j]][i, ] <- sum((betahat != 0) * (simdt$beta != 0))
        
        # Number of false positives
        fpos[[j]][i,] <- nzs[[j]][i,]-tpos[[j]][i,]
        
        # Number of false negatives
        fneg[[j]][i, ] <- sum((betahat == 0) * (simdt$beta != 0))
        
        # Number of true negatives
        tneg[[j]][i, ] <- sum((betahat == 0) * (simdt$beta == 0))
        
        # False positive rate: FPR = FP/(FP + TN)
        fpr[[j]][i, ] <-  fpos[[j]][i, ] / (fpos[[j]][i, ] + tneg[[j]][i, ])
        
        # False negative rate: FNR = FN/(FN + TP)
        fnr[[j]][i, ] <- fneg[[j]][i, ] / (fneg[[j]][i, ] + tpos[[j]][i, ])
        
        # F1 score
        F1[[j]][i, ] <- 2 * tpos[[j]][i, ] / (2 * tpos[[j]][i, ] + 
                                                fpos[[j]][i, ] + fneg[[j]][i, ])
        
        # Matthew correlation coefficient: mltools::mcc(TP =tpos[[j]][i, ], FP=, TN = , FN = )
        numerator <- (tpos[[j]][i, ] * tneg[[j]][i, ]) - (fpos[[j]][i, ] * fneg[[j]][i, ])
        denom <- (tpos[[j]][i, ] + fpos[[j]][i, ]) * (tpos[[j]][i, ] + fneg[[j]][i, ]) * (tneg[[j]][i, ] + fpos[[j]][i, ]) * (tneg[[j]][i, ] + fneg[[j]][i, ])
        # If denom is zero which can happen when TN = 0 OR FN = 0 especially in full model. mcc will be approximately 0 see Chicco and Jurman 2020 
        denominator <- ifelse(denom == 0, 1, denom)
        
        MCC[[j]][i, ] <- numerator/sqrt(denominator)
        
        # Error rate
        error.rate[[j]][i, ] <- (fpos[[j]][i, ] + fneg[[j]][i, ]) / (fpos[[j]][i, ] +
                                                                       fneg[[j]][i, ] + tpos[[j]][i, ] + tneg[[j]][i, ])
        # optimism
        opt[[j]][i, ] <- (err.test[[j]][i, ] - err.train[[j]][i, ]) / err.train[[j]][i, ]
        
        # Relative test error (to the bayes error) (RTE)
        rel.test.err[[j]][i, ] <- err.test[[j]][i, ] / simdt$sigma ^ 2
        
        # Relative Risk (to the null risk) (RR)
        rel.risk[[j]][i, ] <- risk[[j]][i, ] / risk.null[i]

        # estimated betas (without intercept)
        betas[[j]][i, ] <- betahat
        
        # estimated shrinkage factors compared to the least-squares
        shrinkage[[j]][i, ] <- betahat / beta_full
        
        # save tuning parameters
        lambda[[j]][i, ] = reg.obj$lambda
        gamma[[j]][i, ] = reg.obj$gamma
        
      }, error = function(err) {
        if (verbose) {
          cat(paste("Oops! Something went wrong, see error message",
                    "below; recording all metrics here as NAs ...\n"))
          cat("    ***** Error message *****\n")
          cat(sprintf("    %s\n",err$message))
          cat("    *** End error message ***\n")
        }
        # N.B. No need to do anything, the metrics are already filled with NAs
      })
    }
    # Save intermediate results?
    if(!is.null(file) && file.rep > 0 && i %% file.rep == 0) {
      saveRDS(enlistx(err.train,err.val,err.test,err.null,prop,risk,risk.null,
                     nzs,fpos,tpos,fneg,tneg,F1,MCC,error.rate,rel.test.err, opt,
                     sigma,runtime,lambda,gamma,betas,shrinkage,beta.true),
              file = file)
    }
  }
  # Save results now (in case of an error that might occur below)
  out = enlistx(err.train,err.val,err.test,err.null,prop,risk,risk.null,
               nzs,fpos,tpos,fneg,tneg,F1,MCC,fnr,fpr,error.rate,rel.test.err, 
               rel.risk, opt,sigma,lambda,gamma,runtime, betas,shrinkage,
               beta.true)
  if (!is.null(file)) saveRDS(out, file)
  
  out = enlistx(err.train,err.val,err.test,err.null,prop,risk,risk.null,
               nzs,fpos,tpos,fneg,tneg,F1,MCC,fnr,fpr,error.rate,rel.test.err,
               rel.risk, opt,sigma,runtime, betas,shrinkage,lambda,gamma)
  # Save final results
  out = c(out,list(n = n, corrtype = corrtype,betatype = betatype,snr = snr,
                   betatrue = beta.true,n_methods = N,reg.names = reg.names,
                   nrep = nrep, nvars = nvar, call = this.call))
  class(out) = "sim_master"
  if (!is.null(file)) { saveRDS(out, file); invisible(out) }
  else return(out)
}

#' @title Generates correlation structure
#' @description
#' Generates four types of correlation structure as described in Kipruto and Sauerbrei (2022).
#' In addition, it generates block correlation matrix. Block correlation matrix uses the
#' true regression coefficients to determine the indices of signal and noise variables 
#' to create the two groups.
#'@param nvar The number of variables. Default is 15
#'@param type correlation types. `c1` = Houwelingen and Sauerbrei (2013),
#' `c2` = AR1(0.3), `c3` = AR1(0.8), `c4` = Educational Body fat and `c5` = block
#' correlation matrix.
#'@param withinSignalCorr,withinNoiseCorr,betweenSignalNoiseCorr Specifies within signal correlation,
#'within noise correlation and between signal and noise correlation, respectively. 
#'@param betatype Specifies the true regression coefficients. See `beta_type()` function for
#'details.
#' @references Kipruto, E. and Sauerbrei, W. (2022). Comparison of variable 
#' selection procedures and investigation of the role of shrinkage in linear
#' regression-protocol of a simulation study in low-dimensional data. Plos one, 
#' 17(10), p.e0271240.
#'@export 
corr_type <- function(nvar = 15,
                      type = c("c1", "c2", "c3", "c4", "c5"),
                      withinSignalCorr = 0.8, withinNoiseCorr = 0.8,
                      betweenSignalNoiseCorr = 0, betatype = c("a","b","c","d")) {
  
  # match arguments
  type <- match.arg(type)
  betatype <- match.arg(betatype)
  
  if (type=="c1"){
    corrx <- matrix(0, nrow = nvar, ncol = nvar)
    corrx[c(1, 6), c(6, 1)] <- 0.5
    corrx[c(7, 5), c(5, 7)] <- 0.3
    corrx[c(4, 9), c(9, 4)] <- 0.5
    corrx[c(12, 3), c(3, 12)] <- 0.5
    corrx[c(14, 7), c(7, 14)] <- 0.5
    corrx[c(11, 6), c(6, 11)] <- 0.7
    corrx[c(13, 5), c(5, 13)] <- -0.7
    corrx[c(10, 8), c(8, 10)] <- 0.7
    diag(corrx)<- 1
    return(corrx)
  } else if (type=="c2"){
    inds <- 1:nvar
    corrx <- 0.3^abs(outer(inds, inds, "-"))
    return(corrx)
  } else if (type=="c3"){
    inds <- 1:nvar
    corrx<- 0.8^abs(outer(inds, inds, "-"))
    return(corrx)
  } else if (type=="c4"){
    # Initialize the covariance matrix
    corrx <- matrix(0, nrow = nvar, ncol = nvar)
    # Define the correlation matrix for the first 13 variables
    bfat <- matrix(c(1.00,  -0.01,  -0.23, 0.12,  0.17, 0.22, -0.07, -0.20, 0.01, -0.13,  -0.04,   -0.07,  0.22,
                     -0.01,   1.00,  0.52, 0.80,  0.90, 0.87,  0.93,  0.84, 0.83,  0.70,   0.78,    0.75,  0.70,
                     -0.23,   0.52,   1.00, 0.32,  0.26, 0.23,  0.43,  0.34, 0.51,  0.46,   0.31,    0.34,  0.39,
                     0.12,   0.80,   0.32, 1.00,  0.78, 0.74,  0.71,  0.65, 0.64,  0.52,   0.69,    0.71,  0.72,
                     0.17,   0.90,   0.26, 0.78,  1.00, 0.90,  0.81,  0.72, 0.71,  0.56,   0.74,    0.68,  0.65,
                     0.22,   0.87,   0.23, 0.74,  0.90, 1.00,  0.85,  0.73, 0.72,  0.51,   0.67,    0.59,  0.59,
                     -0.07,   0.93,   0.43, 0.71,  0.81, 0.85,  1.00,  0.88, 0.80,  0.63,   0.74,    0.68,  0.59,
                     -0.20,   0.84,   0.34, 0.65,  0.72, 0.73,  0.88,  1.00, 0.77,  0.61,   0.74,    0.67,  0.50,
                     0.01,   0.83,   0.51, 0.64,  0.71, 0.72,  0.80,  0.77, 1.00,  0.73,   0.63,    0.63,  0.65,
                     -0.13,   0.70,   0.46, 0.52,  0.56, 0.51,  0.63,  0.61, 0.73,  1.00,   0.53,    0.56,  0.65,
                     -0.04,   0.78,   0.31, 0.69,  0.74, 0.67,  0.74,  0.74, 0.63,  0.53,   1.00,    0.76,  0.60,
                     -0.07,   0.75,   0.34, 0.71,  0.68, 0.59,  0.68,  0.67, 0.63,  0.56,   0.76,    1.00,  0.64,
                     0.22,   0.70,   0.39, 0.72,  0.65, 0.59,  0.59,  0.50, 0.65,  0.65,   0.60,    0.64,  1.00),
                   ncol = 13, nrow = 13)
    
    # Assign the correlation matrix to the corresponding section in the final matrix
    corrx[1:13, 1:13] <- bfat
    
    # Generate the remaining uncorrelated variables
    if (nvar > 13) {
      diag(corrx)[14:nvar] <- 1
    }
  } else {
    corrx <- matrix(0, nrow = nvar, ncol = nvar)
    
    # Use betatype to identify the position of signal and noise variables
    truebeta <- beta_type(nvar = nvar, type = betatype)
    signal_indices <- which(truebeta!=0)
    noise_indices <- setdiff(1:nvar, signal_indices)
    
    # Fill within-group correlations. We assume the first 7 variables are signal and the
    # rest are noise
    corrx[signal_indices, signal_indices] <- withinSignalCorr
    corrx[noise_indices, noise_indices] <- withinNoiseCorr
    
    # Fill between-group correlations
    corrx[signal_indices, noise_indices] <- betweenSignalNoiseCorr
    corrx[noise_indices, signal_indices] <- betweenSignalNoiseCorr
  }
  return(corrx)
  
}

#' @title Generates true regression coefficients
#' @description
#' Generates four types of regression coefficients
#' @param nvar The number of variables. Default is 15
#' @param type The type of regression coefficients. See Table 1 of the reference
#' @references Kipruto, E. and Sauerbrei, W. (2022). Comparison of variable 
#' selection procedures and investigation of the role of shrinkage in linear
#' regression-protocol of a simulation study in low-dimensional data. Plos one, 
#' 17(10), p.e0271240.
#'@export 
beta_type <- function(nvar = 15,
                      type = c("a", "b", "c", "d")) {
  
  # Match arguments
  type <- match.arg(type)
  
  if (type=="a"){
    out <- c(1.5, 0, 1, 0, 1, 0, 0.5, 0, 0.5, 0, 0.5, 0,-0.5, 0, 0)
    beta <-replace(rep.int(0, times = nvar), list = 1:15, values =out)
    return(beta)
  } else if (type=="b"){
    out <- c(1.5, 0, 0.5, 0, 0.5, 0, 0.25, 0, 0.25, 0, 0.25, 0,-0.25, 0, 0)
    beta <-replace(rep.int(0, times = nvar), list = 1:15, values =out)
    return(beta)
  } else if (type=="c"){
    index <- floor(seq(1,nvar,nvar/7))
    beta <-replace(rep(0, times = nvar), list = index, values =1)
    return(beta)
  } else {
    beta <-replace(rep(0, times = nvar), list = seq(1,7), values =1)
    return(beta)
  }
}


#' Enlist function
#' @param ... parameters of sim_master
#' @export
enlistx <- function (...)
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}

#' @title Standard Deviation of a Vector
#' @description
#' Calculates the standard deviation using 'n' instead of '(n-1)' in the denominator,
#' aligning with the conventions of the glmnet package.
#' @param x A vector of numeric values. Missing data or `NA` will be removed.
#' @export
sd_x <- function(x) {
  x <- x[!is.na(x)]
  meanx <- mean(x)
  n <- length(x)
  sqrt(sum((x - meanx) ^ 2) / n)
}

#' @title Estimate residual variance from OLS model
#' @description
#' Calculates the residual variance (mean squared error) from a 
#' least-squares model using a standardized matrix of predictors and a centered
#' response variable.
#' @param x A standardized matrix of predictors.
#' @param y A centered response variable.
#' @return An estimate of the residual variance (`sigma2`) and regression 
#' estimates (`coefficients`) without intercept.
#' @export
residual_variance <- function(x, y) {
  
  # Number of variables and observations
  x <- as.matrix(x)
  n <- dim(x)[1L]
  
  # Fit linear model
  fit <- .lm.fit(y = y, x = x)
  
  # Estimate model residual variance. Intercept assumed to be zero due to centering of y and 
  # standardization of x
  sigma2 <- sum(fit$residuals ^ 2) / (n - fit$rank)
  
  # Parameters of interest: sigma2 and regression coefficients
  model_parameters <- list(sigma2 = sigma2,
                           coefficients = fit$coefficients)
  return(model_parameters)
}
