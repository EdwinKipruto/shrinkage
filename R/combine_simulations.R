#' Combine simulation output
#'
#' Combine simulation output stored in an object of class 'sim_master' (produced
#' by \code{\link{sim_master}}) into a dataframe for further processing. This
#' function is not suitable for regression estimates or shrinkage factors.
#' Please refer to the 'combine_sim_coef' function for those cases.
#'
#' @param path  the directory where the rds files are stored for each simulation
#' scenario. E.g "C:/Users/ek/Desktop/varselection/"
#' @param metric  the metric of interest. Available metrics include:
#' @param add_rep Specifies whether to include simulation replication in the output. Important for
#'  calculating similarities. Set it to TRUE if `compute_similarity_sim()` function will be used.
#' \itemize{
#'   \item `rr` - relative risk which measures the risk or model error relative to the null model (model without predictors).
#'   \item `rte` - relative test error which measures the expected test error relative to the Bayes error rate.
#'   \item `trainerror` - error on the training data.
#'   \item `valerror` - Error on the validation data.
#'   \item `testerror` - error on the test data.
#'   \item `nullerror` - error for a model without predictors.
#'   \item `prop` - proportion of variance unexplained by the model.
#'   \item `risk` - model error.
#'   \item `nullrisk` - error in the null model.
#'   \item `nonzero` - Number of Nonzero Coefficients.
#'   \item `fpos` - false Positives.
#'   \item `tpos` - true Positives.
#'   \item `fneg` - false Negatives.
#'   \item `tneg` - True Negatives.
#'   \item `F1` - F-score which measures the accuracy of the support recovery.
#'   \item `MCC` - Matthews correlation coefficient which measures the accuracy of the support recovery.
#'   \item `fpr` - False Positive Rate which measures the proportion of actual negatives that are predicted as positives.
#'   \item `fnr` - False Negative Rate which measures the proportion of actual positives that are predicted as negatives.
#'   \item `errorrate` - Error rate.
#'   \item `optimism` - Optimism of training error is computed as (test error - train error)/train error.
#'   \item `runtime` - Time used in simulation runs.
#'   \item `lambda` - Tuning Parameter for Regression Models.
#'   \item `gamma` - Tuning Parameter for Regression Model.
#'}
#'
#' @return A dataframe containing the combined simulation output for the specified metric.
#'
#' @export

combine_sim_metrics <- function(path,
                                metric = c("rr", "rte", "trainerror",
                                           "valerror", "testerror",
                                           "nullerror", "prop", "risk",
                                           "nullrisk", "nonzero", "fpos",
                                           "tpos", "fneg", "tneg", "F1","MCC",
                                           "fpr", "fnr", "errorrate",
                                           "optimism", "lambda", "gamma",
                                           "runtime"), add_rep = FALSE){

  # check whether the path to the folders where rds files are stored is valid
  if (!dir.exists(path))
    stop("The path is invalid", call. = FALSE)

  # match arguments
  metric <- match.arg(metric)

  # The path of each RDS files in the directory
  files <- paste0(path, list.files(path, pattern = "\\.rds$"))

  # Read RDS files and store data in a list: length of list = number of scenarios
  data_list <- lapply(files, function(v) readRDS(v))

  # total number of scenarios
  nscenarios <- length(data_list)

  # number of replications. we assume they are fixed for all scenarios
  nrep <- data_list[[1]]$nrep

  if (is.null(nrep)) {
    stop(sprintf("Number of replications (%d) not found in the RDS files", nrep))
  }

  # simulation parameters for each scenario for each replication
  repn <- 1:nrep
  if(add_rep){
    sim_parms <- lapply(data_list, function(v) data.frame(rep = repn,
                                                          betatype = rep(v$betatype,nrep),
                                                          corrtype = rep(v$corrtype,nrep),
                                                          snr = rep(v$snr,nrep),
                                                          n = rep(v$n,nrep),
                                                          sigma = v$sigma,
                                                          nvars = rep(v$nvars,nrep))

    )
  }else{
    sim_parms <- lapply(data_list, function(v) data.frame(betatype = rep(v$betatype,nrep),
                                                          corrtype = rep(v$corrtype,nrep),
                                                          snr = rep(v$snr,nrep),
                                                          n = rep(v$n,nrep),
                                                          sigma = v$sigma,
                                                          nvars = rep(v$nvars,nrep))
    )

  }

  # extract the metrics of interest: each column is a metric for a method
    sim_metrics <- lapply(data_list, function(v) do.call(cbind,
                  switch (metric,
                          "trainerror" = v$err.train,
                          "valerror" = v$err.val,
                          "testerror" = v$err.test,
                          "nullerror" = v$err.null,
                          "prop" = v$prop,
                          "risk" = v$risk,
                          "nullrisk" = v$risk.null,
                          "nonzero" = v$nzs,
                          "fpos" = v$fpos,
                          "tpos" = v$tpos,
                          "fneg" = v$fneg,
                          "tneg" = v$tneg,
                          "F1" = v$F1,
                          "MCC" = v$MCC,
                          "fnr" = v$fnr,
                          "fpr" = v$fpr,
                          "errorrate" = v$error.rate,
                          "rte" = v$rel.test.err,
                          "rr" = v$rel.risk,
                          "optimism" = v$opt,
                          "lambda" = v$lambda,
                          "gamma" = v$gamma,
                          "runtime" = v$runtime
                          )
                  )
                  )

    # assign names to the column of the matrices
    # extract names of regression methods
    cnames <-lapply(data_list,function(v) v$reg.names)

    # assign names to metrics of interest
    sim_metrics <- Map(function(datax, xnames) {
           colnames(datax) <- xnames
           datax
       }, sim_metrics, cnames)


    # combine simulation parameters and training errors for each scenario
    bind <- Map(cbind, sim_parms, sim_metrics)

    # combine list into large datatable.
    fmat <- data.table::rbindlist(bind, fill = TRUE)

    # add scenario
    #fmat <- cbind(scenario = rep(1:nscenarios, each = nrep), fmat)
    fmat <- data.table(scenario = rep(1:nscenarios, each = nrep), fmat)

    # attach an attribute: important for plotting
    attr(fmat, "metric") <- metric
    attr(fmat, "add_rep") <- add_rep

    # to use dplyr::filter we need a data.frame as class
    class(fmat) <- c("data.table", "data.frame", "combinesim")
    return(fmat)
}

#' Combines regression coefficients or shrinkage factors
#'
#' #The function combines the output stored in an object of class sim_master
#' (produced by \code{\link{sim_master}}) into a dataframe. In cases where
#' different numbers of variables are used, the function combines the inputs
#' into one large data frame. If there are missing columns in the smaller data
#' frames, the function fills those columns with NAs in the combined data frame.
#'
#' @param path The directory where the RDS files for each simulation scenario
#' are stored. Example: "C:/Users/Desktop/varselection/"
#' @param type Specifies whether to return "betas" or "shrinkage". Default is "betas",
#' which returns a dataframe of regression estimates for each method and
#' replication of a scenario.
#' @param freq A logical value that indicates whether the regression estimates of
#' variables should be converted to 0s and 1s. This is useful if the inclusion
#' frequency of variables is required. If \code{freq = TRUE} and
#' \code{plot_coefs()} is used, the inclusion frequencies will be plotted.
#' @param add_rep Specifies whether to include simulation replication in the output. Important for
#'  calculating similarities. Set it to TRUE if `compute_similarity_sim()` function will be used.
#' @return A dataframe for each scenario without aggregation. This means that
#' each scenario is repeated for each replication. If the number of replications
#' is 5, then each scenario will be replicated 5 times for each regression method.
#' @export
combine_sim_coef <- function(path, type = c("betas", "shrinkage"), add_rep = FALSE, freq = FALSE) {


  # check whether the path to the folders where RDS files are stored is valid
  if(!dir.exists(path)) {
    stop("The path is invalid", call. = FALSE)
  }

  # match argument
  type <- match.arg(type)

  # The path of each RDS files in the directory
  files <- paste0(path,list.files(path, pattern = "\\.rds$"))

  # Read RDS files and store data in a list: length of list = number of scenarios
  data_list <- lapply(files, function(v) readRDS(v))

  # number of replications and methods used in the simulation. we assume the
  # number of replications are fixed for all scenarios
  nrep <- unlist(lapply(data_list, function(v) v$nrep), use.names = FALSE)

  if (length(unique(nrep)) != 1)
     warning("Different numbers of repetitions were used. The program assumes that the number of repetitions should be the same in all simulations. Only the first repetition is used.
            ")
  nrep <- nrep[1]

  #N <- data_list[[1]]$n_methods
  N <- lapply(data_list, function(v) v$n_methods)

  # assert that nrep and method must be in the files
  if (is.null(nrep)) {
    stop("number of replications (nrep) not found in the files", call. = FALSE)
  }

  if (any(is.null(N))) {
    stop("At least one element in the number of methods (n_methods) is null in the files", call. = FALSE)
  }
  # total number of scenarios
  nscenarios <- length(data_list)

  # simulation parameters for each scenario for each replication
  repn <- 1:nrep
  if(add_rep){
  sim_parms <- lapply(data_list, function(v) data.frame(rep = repn,
                                                        betatype = rep(v$betatype,nrep),
                                                        corrtype = rep(v$corrtype,nrep),
                                                        snr = rep(v$snr,nrep),
                                                        n = rep(v$n,nrep),
                                                        sigma = v$sigma,
                                                        nvars = rep(v$nvars,nrep))

  )
  }else{
    sim_parms <- lapply(data_list, function(v) data.frame(betatype = rep(v$betatype,nrep),
                                                          corrtype = rep(v$corrtype,nrep),
                                                          snr = rep(v$snr,nrep),
                                                          n = rep(v$n,nrep),
                                                          sigma = v$sigma,
                                                          nvars = rep(v$nvars,nrep))
    )

  }

  # regression estimates and shrinkage factors
  sim_metrics <- lapply(data_list, function(v) if (type=="betas") {
                                                     v$betas} else {
                                                               v$shrinkage
                                                       }
                          )

  # Assign column names to each matrix of betas/shrinkage in the nested list
  sim_metrics <- lapply(sim_metrics, function(v) {
    lapply(v, assign_column_names)
  }
  )

    # Column-bind corresponding lists for simulation parameters ('sim_parms)
    # and estimated regression coefficients or shrinkage factors ('sim_metrics')
    flist <- Map(function(x, y) {
                                 lapply(x, function(v) cbind(y,v))
                                 }, x=sim_metrics, y=sim_parms
    )
    # rowbind betas/shrinkage for each method in each scenario
    scenario_specific_parms <-lapply(flist,
                                     function(v)data.table::rbindlist(v, fill=TRUE,
                                                                      idcol = "method"))
    # row bind all scenarios together and return a large data.frame. If the number
    # of variables are of different length fill with NAs the missing columns
    #all_scenarios_parms=data.table::rbindlist(scenario_specific_parms, fill=T, idcol=T)
    dfx = data.table::rbindlist(scenario_specific_parms, fill=T, idcol=F)

    # add scenario for easy aggregation
    N <- unlist(N, use.names = TRUE)
   # dfx <- cbind(scenario = rep(1:nscenarios, times = nrep*N), dfx)
    dfx[, scenario := rep(1:nscenarios, times = nrep * N)]

    # order columns for easy processing
    xnames <- c("scenario", "betatype", "corrtype",  "snr", "n", "sigma", "nvars", "method")
   if(add_rep){
    xnames <- c("scenario","rep", "betatype", "corrtype",  "snr", "n", "sigma", "nvars", "method")
   }
    index <- which(!colnames(dfx)%in%xnames)
    data.table::setcolorder(dfx, c(xnames, colnames(dfx)[index]))

    # Place the column of "methods" after the "nvars" column for easier processing later.
    #dfx <- dplyr::relocate(dfx, method, .after = nvars)

    # convert data to 0s and 1s for calculation of inclusion frequency
    if (freq) {
     cols_to_update <- if (!add_rep) 9:ncol(dfx) else 10:ncol(dfx)
     dfx[, (cols_to_update) := lapply(.SD, function(x) ifelse(is.na(x), NA, ifelse(x == 0, 0, 1))), .SDcols = cols_to_update]
    }
    # add attribute for plotting purpose in plot_coefs()
    attr(dfx, "metric") <- ifelse(freq, "inclusion", type)
    attr(dfx, "add_rep") <- add_rep
    # Assign a class to the object. Use "data.frame" for important operations,
    # while use "combinecoef" for easy identification in plotting.
    class(dfx) = c(class(dfx),"combinecoef")
    return(dfx)
}


#' Helper function to assign column names to a matrix
#'
#' @param x the matrix to which column names will be assigned
#' @export
assign_column_names <- function(x) {
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  colnames(x) <- paste0("x", 1:dim(x)[2L]) #sprintf("x%d", 1:dim(x)[2L])
  return(x)
}

#' Extract Simulation Scenario Data from RDS Files
#'
#' This function reads multiple RDS files containing simulation scenario data,
#' extracts relevant information, and combines it into a single data frame.
#'
#' @param files_paths A character vector of file paths to the RDS files.
#'
#' @return A data frame containing the extracted information from all the RDS files.
#' Each row corresponds to a simulation scenario and includes the following columns:
#' \describe{
#'   \item{n}{The sample size of the simulation scenario.}
#'   \item{snr}{The signal-to-noise ratio of the simulation scenario.}
#'   \item{betatype}{The type of beta used in the simulation scenario.}
#'   \item{corrtype}{The type of correlation used in the simulation scenario.}
#'   \item{nvar}{The number of variables used in the simulation scenario.}
#' }
#' @examples
#' \dontrun{
#' # Assuming you have RDS files are stored in the path
#' path <- "Q:/x/postestimation/simdata/final/"
#' sim_data <- extract_sim_scenario(paths)
#' print(sim_data)
#' }
#' @export
extract_sim_scenario <- function(files_paths) {
  path <- list.files(path = files_paths, pattern = "\\.rds$", full.names = TRUE)
  datalist <- lapply(path, function(file) {
    data_scenaria <- readRDS(file)
    data.frame(
      n = data_scenaria$n,
      snr = data_scenaria$snr,
      betatype = data_scenaria$betatype,
      corrtype = data_scenaria$corrtype,
      nvar = length(data_scenaria$betatrue)
    )
  })

  return(do.call(rbind, datalist))
}

