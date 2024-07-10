#' Plot the mean model error,relative test error, etc over several simulation settings.
#'
#'@param x The output of the `combine_sim_metrics()` function.
#'@param betatypes the beta type used in the simulation. Default is `betatypes = "a"`
#'@param corrtypes the correlation type type used in the simulation.
#'Default is `corrtypes = "a"`.
#'@param samplesize a vector of sample size used in the simulation. Default is NULL
#'and all the sample sizes used are plotted.
#'@param snrx a vector of signal to noise ratio used in the simulation. Default is NULL
#'and all the snr used are plotted.
#'@param nvar a numeric value for the number of variables used in the simulation.
#'Default is 15 variables
#'@param methods the name(s) of the method used in the simulation. Default is
#'null and all methods are plotted.
#'@param xaxis parameter to plot in the x axis. Default is `xaxis = snr`
#'@param plotype the plot to generate. Default is `plotype = "means"` which plots
#'the average values of a given metric. If distributions are required,
#'set `plotype = "distribution"` to generate boxplots.
#'@param legend.position the position of the legend in the plot. Default is to
#'place it on the right side of the plot
#'@param row,col One of `betatype` or `corrtype`,indicating the variables for
#'the rows and columns for plotting facet grid; note that row and col must be
#'different. Default is `row = betatype` and `col = corrtype`, so that the
#'plotting grid displays the given metric versus the SNR or sample size,
#'in a plotting grid with the coefficient types across the rows, and the
#'correlation levels across the columns.
#'@param linetype Specifies whether points with corresponding lines should be
#'plotted for the averages. The default is FALSE, which means only averages are
#'plotted. Please note that this option is only relevant when `plotype = "means"`.
#'@param add.se Specifies whether to include one standard error (SE) in
#'addition to the average values of the given metric. The default is FALSE. This
#'option is only relevant when `plotype = "means"`.
#'@param logx specifies whether to log-transform the mean when generating the plot
#'@param lwd,main,ylim,ylab,xlab,facet graphical parameters. see `[ggplot()]` and
#' `[plot]` for details
#'@param boxwidth boxplot width
#'@param points_size size of points in the plot
#'@param points_shape shape of points in the plot
#'@param dodgewidth determines the amount of separation between the grouped elements in
#'boxplot
#' @param facet_size Text size of the facet.
#' @param facet_color Text color of the facet.
#' @param recode_methods Specifies whether to recode the names of the methods based on the
#' user-supplied 'new_methods_names'
#' @param new_methods_names List of new names for the method which is used by recode() function
#' in dplyr. For example, if the old names are c("lasso(cv)", "lasso(aic)) we can change the names
#' to just cv and aic using this approach: list(`lasso(cv)` = "cv", `lasso(aic)` = "aic")
#' @param legend_title The name of the legend title. Suitable when new methods_name is used description
#' @param add.xlab,add.ylab specifies whether to add x and y labels
#'@import data.table
#'@import ggplot2
#'@export
plot_metric <- function(x,
                  betatypes=c("a", "b", "c", "d"),
                  corrtypes=c("c1", "c2", "c3", "c4", "c5"),
                  samplesize = NULL,
                  snrx = NULL,
                  nvar = 15,
                  methods = NULL,
                  xaxis = c("snr", "samplesize"),
                  plotype = c("means", "median", "distribution"),
                  legend.position = "right",
                  row=c("betatype", "corrtype"),
                  col=c("corrtype", "betatype"),
                  linetype = FALSE,
                  add.se = FALSE,
                  lwd = 0.7,
                  points_shape = 16,
                  points_size = 2,
                  dodgewidth = 0.5,
                  main = NULL,
                  ylim = NULL,
                  ylab = NULL,
                  xlab = NULL,
                  boxwidth = 0.5,
                  logx = FALSE,
                  facet = FALSE,
                  facet_size = 13,
                  facet_color = "black",
                  recode_methods = FALSE,
                  new_methods_names = NULL,
                  legend_title = NULL,
                  add.xlab = TRUE,
                  add.ylab = TRUE
                  ){
  # Check for ggplot2 package
  if (!requireNamespace("ggplot2",quietly=TRUE)) {
    stop("Package ggplot2 not installed. Please install it!")
  }

 # calculate the metric of interest if x is path where the rds files are stored
 #x <- combine_simulations(path, metric = type)

  # assert that x must be of class "combinesim"
  if (!inherits(x, "combinesim")) {
    stop("x must be an output from combine_sim_metrics() function",call. = FALSE )
  }

  # match arguments
  corrtypes <- match.arg(corrtypes)
  betatypes <- match.arg(betatypes)
  #legend.position <- match.arg(legend.position)

  col<-match.arg(col)
  row<-match.arg(row)
  plotype <- match.arg(plotype)
  xaxis<-match.arg(xaxis)

  # Assert that row and col must be different
  if (row==col) {
    stop("row and col must be different", call. = FALSE)
  }

  # metric used in combine_sim_metrics()
  type <- attr(x, "metric")

  if(!corrtypes%in%unique(x$corrtype)){
    stop(sprintf("The corrtypes = %s not used in the simulation.", corrtypes), call. = FALSE)

  }

  if(!betatypes%in%unique(x$betatype)){
    stop(sprintf("The betatypes = %s not used in the simulation.", betatypes), call. = FALSE)
  }

  if(!is.null(samplesize) && !all(samplesize%in%unique(x$n))){
    stop(sprintf("Not all sample sizes supplied were used in the simulation. This applies to sample size:%s",
                 paste0(setdiff(samplesize, unique(x$n)), collapse = ", ")), call. = FALSE)
  }

  if(!is.null(snrx) && !all(snrx%in%unique(x$snr))){
    stop(sprintf("Not all SNR supplied were used in the simulation. This applies to SNR:%s",
            paste0(setdiff(snrx, unique(x$snr)), collapse = ", ")), call. = FALSE)
  }

  # Set y-label
  if (is.null(ylab)& add.ylab) {
  ylab  = switch(type,
                 "trainerror" = "Training Error",
                 "valerror"   = "Validation Error", # what is the difference btwn val and test?
                 "testerror"  = "Test Error",
                 "nullerror"  = "Null Error",
                 "prop"       = "Proportion Variance",
                 "risk"       = "Model Risk",
                 "nullrisk"   = "Null Risk",
                 "rte"        = "Relative test error(to Bayes)",
                 "errorrate"  = "Error rate",
                 "rr"         = "Relative Risk",
                 "optimism"   = "Model Optimism",
                 "nonzero"    = "Number of Nonzeros",
                 "fpos"       = "False Positive Rate",
                 "tpos"       = "True Positive Rate",
                 "fneg"       = "False Negative Rate",
                 "tneg"       = "TNR",
                 "fpr"        = "FPR",
                 "fnr"        = "FNR",
                 "F1"         = "F1 score",
                 "MCC"         = "Matthew Correlation",
                 "lambda"     = "Lambda",
                 "gamma"      = "Gamma"
)
}
  # all methods used in the analysis: column no 8+. If methods for plotting are
  #not provided we plot all methods
  methodx <- colnames(x)[-c(1:7)]

  if (!is.null(methods)) {
          invalid_methods <- methods[!methods %in% methodx]
     if (length(invalid_methods) > 0){
       stop("The methods ", invalid_methods, " are invalid", call. = FALSE)
     }
    # use user methods
    methodx <- methods
  }

  # Deal with signal to noise ratio
  snr_used <- unique(x$snr)

  if (!is.null(snrx)) {
    snr_index <- which(!snrx%in%snr_used)
    if (length(snr_index)!=0) {
      stop("The snr supplied was not used in the simulation. \n",
           sprintf("i This includes snr of: %s.",
                   paste0(snrx[snr_index], collapse = ", ")))
    }
    # snr to be used for plotting
    snr_used <- unique(snrx)
  }

  # Deal with sample size
  n_used <- unique(x$n)
  if (!is.null(samplesize)) {
       ind <- which(!samplesize%in%n_used)
        if (length(ind)!=0){
            stop("The sample supplied was not used in the simulation. \n",
               sprintf("i This includes sample size of: %s.",
                                paste0(samplesize[ind], collapse = ", ")))
        }
  n_used <- unique(samplesize)
  }

  # Deal number of covariates used in the model
  if (length(nvar)!=1) {
    stop("! nvar must be a single number specifying the number of variables used.\n",
         sprintf("i The length of nvar is %d.", length(nvar)))
  }

  nvar_used <- unique(x$nvars)
  nvar_index <- which(!nvar%in%nvar_used)

  if (length(nvar_index)!=0){
    stop("! nvar = ", nvar, " was not used in the simulation.", call. = FALSE)
  }

  # remove missing and subset the data based on simulation parameters and convert from wide to long
  #x <- data.table::as.data.table(x)
  # select methods before omitting NA: Added 17.01.2024
  select_columns <- c(colnames(x)[c(1:7)], methodx)
  x <- x[, select_columns, with = FALSE]
  x <- x[complete.cases(x)]

  if(nrow(x) == 0){
    stop("x is empty after removing NA", call. = FALSE)
  }

  me1<- x[betatype==betatypes & corrtype==corrtypes &snr%in%snr_used & n%in%n_used & nvars==nvar]
  # convert the data to data.table to be used by melt() function for long format
  long_data <- data.table::melt(me1,
                                id.vars = c("scenario", "betatype","corrtype",
                                            "snr", "n","sigma", "nvars"),
                                variable.name = "method",
                                value.name = "metric")
  # subset the methods
  long_data <- long_data[method%in%methodx]
  # calculate the mean within a group with corresponding standard error and add
  # to the long_data by replicating the values within a group
  dfx <- long_data[, c(.SD,
                       median_value = smedian(metric)[[1]],
                       median_lower = smedian(metric)[[2]],
                       median_upper = smedian(metric)[[3]],
                       mean_value = if(logx){log(mean(metric))}else{mean(metric)},
                       se = std_mean(metric)),
                   by = .(scenario,method)]

  # Recode methods if necessary like lasso(cv) = cv, lasso(aic) = aic etc
  if (recode_methods){
    if (is.null(new_methods_names)){
      stop("Supply new methods for recoding through new_methods_names argument", call. = FALSE)

    } else {
      if (!any(names(new_methods_names)%in%unique(dfx$method))){
        stop(sprintf("All new_methods_names are not part of the methods: %s",unique(dfx$method)), call. = FALSE)
      }
    }
    dfx <- dplyr::mutate(dfx, method = dplyr::recode(dfx$method, !!!new_methods_names))
  }

  # This part is relevant when methods are recoded
  if (is.null(legend_title)){
    legend_title <- "method"
  }
  # change the name of the variable which can also change the legend title
  setnames(dfx, "method", legend_title)


  # rename betatype and corrtypes for plotting
  dfx[betatype == betatypes, betatype := paste0("Beta-type = ", betatypes)]
  dfx[corrtype == corrtypes, corrtype := paste0("Corr-type = ", corrtypes)]
  # important for plotting
  dfx$method <- factor(dfx$method, levels=unique(dfx$method))



  # Deal with x axis and subset data based on the x axis option
  if (xaxis=="snr"){
    if(is.null(xlab) & add.xlab){
    xlab <- "Signal-to-noise ratio"
    }
    sampx <- n_used
    nx <- length(sampx)
  } else {
    if (is.null(xlab)& add.xlab){
    xlab <- "Sample Size"
    }
    sampx <- snr_used
    nx <- length(sampx)
  }
  # Save plots as list
  plots <- vector("list", nx)

  if(plotype=="means"){
  # Get rid of the metric that has been replaced by the statistic.
  # Important for removing duplicates
  dfx <- subset(dfx, select = -metric)
  dfx <- dfx[!duplicated(dfx), ]

  for (v in 1:nx) {
    if (xaxis == "snr") {
      datax <- dfx[n == sampx[v]]
      mainx <- paste0("n = ", sampx[v])
      plots[[v]] <- ggplot2::ggplot(data = datax,
                                    ggplot2::aes_string(x = "factor(.data$snr)",
                                                        y = ".data$mean_value",
                                                        color = paste0(".data$", legend_title),
                                                        linetype = paste0(".data$", legend_title),
                                                        group = paste0(".data$", legend_title)))  +
        ggplot2::geom_point(position = position_dodge(width = dodgewidth), shape = points_shape, size = points_size) +
        theme_bw() +  ggplot2::theme(legend.position = legend.position) +
        ggplot2::labs(x = xlab, y = ylab)
        # if(type != "nonzero"){
        #   plots[[v]] <- plots[[v]] + ggplot2::scale_y_continuous(labels = scaleFUN)
        # }
      if (add.se) plots[[v]] <- plots[[v]] +
        geom_errorbar(aes(ymin = .data$mean_value - .data$se, ymax = .data$mean_value + .data$se),
                    width = 0.2, position = position_dodge(width = dodgewidth))
      # linetype might not connect if many methods are supplied
      if (linetype) plots[[v]] <- plots[[v]] + ggplot2::geom_line(ggplot2::aes_string(x = "factor(.data$snr)",
                                                                      y = ".data$mean_value",
                                                                      linetype = paste0(".data$", legend_title),
                                                                      color = paste0(".data$", legend_title),
                                                                      group = paste0(".data$", legend_title)),position = position_dodge(width = dodgewidth), lwd = lwd)

    } else {
      datax <- dfx[snr == sampx[v]]
      datax$n <- factor(datax$n)
      mainx <- paste0("snr = ", sampx[v])

      plots[[v]] <- ggplot2::ggplot(data = datax,
                                    ggplot2::aes(x = factor(.data$n),
                                                 y = .data$mean_value,
                                                 color = .data$method)) +
        ggplot2::geom_point(position = position_dodge(width = dodgewidth), shape = points_shape, size = points_size) +
        theme_bw() +  ggplot2::theme(legend.position = legend.position) +
        ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
      if (add.se) plots[[v]] <- plots[[v]] +  geom_errorbar(aes(ymin = .data$mean_value - .data$se, ymax = .data$mean_value + .data$se),
                                                            width = 0.2, position = position_dodge(width = dodgewidth))
      if (linetype) plots[[v]] <- plots[[v]] + ggplot2::geom_line(aes(group = method, linetype = method),position = position_dodge(width = dodgewidth), lwd = lwd)
      #plots[[v]] <-  plots[[v]] + guides(group = guide_legend(title = legend_title))

      }

    if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col))) +
        theme(strip.text = element_text(size = facet_size, color = facet_color))
    if (is.null(main)) plots[[v]] <- plots[[v]] + ggplot2::ggtitle(mainx)
    if (!is.null(ylim)) plots[[v]] <- plots[[v]] + ggplot2::coord_cartesian(ylim = ylim)

  }

  } else if(plotype=="median") {

    dfx <- subset(dfx, select = -metric)
    dfx <- dfx[!duplicated(dfx), ]

    for (v in 1:nx) {
      if (xaxis == "snr") {
        datax <- dfx[n == sampx[v]]
        mainx <- paste0("n = ", sampx[v])
        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = factor(.data$snr),
                                                   y = .data$median_value,
                                                   color = .data$method)) +
          ggplot2::geom_point(position = position_dodge(width = dodgewidth), shape = points_shape, size = points_size) +
          theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
        if (add.se) plots[[v]] <- plots[[v]] +
          geom_errorbar(aes(ymin = .data$median_lower, ymax = .data$median_upper),
                        width = 0.2, position = position_dodge(width = dodgewidth))
        if (linetype) plots[[v]] <- plots[[v]] + ggplot2::geom_line(aes(linetype = method, color = method, group = method),position = position_dodge(width = dodgewidth), lwd = lwd)
        plots[[v]] + guides(color = guide_legend(title = legend_title))

      } else {
        datax <- dfx[snr == sampx[v]]
        datax$n <- factor(datax$n)
        mainx <- paste0("snr = ", sampx[v])

        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = factor(.data$n),
                                                   y = .data$median_value,
                                                   color = .data$method)) +
          ggplot2::geom_point(position = position_dodge(width = dodgewidth), shape = points_shape, size = points_size) +
          theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + guides(color = guide_legend(title = legend_title))
        if (add.se) plots[[v]] <- plots[[v]] +  geom_errorbar(aes(ymin = .data$median_lower, ymax = .data$mean_upper),
                                                              width = 0.2, position = position_dodge(width = dodgewidth))
        if (linetype) plots[[v]] <- plots[[v]] + ggplot2::geom_line(aes(group = method, linetype = method),position = position_dodge(width = dodgewidth), lwd = lwd)
        plots[[v]] + guides(color = guide_legend(title = legend_title))

      }

      if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col))) +
          theme(strip.text = element_text(size = facet_size, color = facet_color))
      if (is.null(main)) plots[[v]] <- plots[[v]] + ggplot2::ggtitle(mainx)
      if (!is.null(ylim)) plots[[v]] <- plots[[v]] + ggplot2::coord_cartesian(ylim = ylim)

    }


  } else {
    # plot distributions
    for (v in 1:nx) {
      if (xaxis == "snr") {
        datax <- dfx[n == sampx[v]]
        datax$snr <- as.factor(datax$snr)
        mainx <- paste0("n = ", sampx[v])
        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = .data$snr,
                                                   y = .data$metric,
                                                   colour = .data$method)) +
          ggplot2::geom_boxplot(width = boxwidth) +theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + guides(color = guide_legend(title = legend_title))
      } else {
        datax <- dfx[snr == sampx[v]]
        datax$n <- as.factor(datax$n)
        mainx <- paste0("snr = ", sampx[v])

        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = .data$n,
                                                   y = .data$metric,
                                                   colour = .data$method)) +
          ggplot2::geom_boxplot(position=position_dodge(width=dodgewidth), width = boxwidth) + theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + guides(color = guide_legend(title = legend_title))
      }

      if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col))) +
          theme(strip.text = element_text(size = facet_size, color = facet_color))
      if (is.null(main)) plots[[v]] <- plots[[v]] + ggplot2::ggtitle(mainx)
      if (!is.null(ylim)) plots[[v]] <- plots[[v]] + ggplot2::coord_cartesian(ylim = ylim)

  }

}
  plots
}

#' Scale Function
#'
#' This function scales numeric values to a string representation with two decimal places.
#'
#' @param x A numeric value.
#' @return A string representation of the input numeric value with two decimal places.
scaleFUN <- function(x) {
  sprintf("%0.3f", x)
}


#' @title Plots regression estimates, shrinkage factors or inclusion frequencies
#' @description
#' This function generates plots for the regression estimates, shrinkage factor,
#' or inclusion frequency of each variable. It provides the option to plot either
#' the means or distributions.
#' @param x The output of the \code{combine_sim_coef()} function. If inclusion
#' frequencies are required, the user must set \code{freq = TRUE} in
#' \code{combine_sim_coef()} and \code{plotype = "means"} in this function.
#' @param truncate.lower,truncate.upper Specifies whether to truncate the
#' lower shrinkage factors to truncate.lower and/or the upper shrinkage factors
#' to truncate.upper. The default is NULL, which means no truncation is conducted.
#' It is important to note that truncation is performed before calculating the
#' averages and standard errors.
#' @param varx The names of variables to plot. The default value is NULL,
#' indicating that all variables will be plotted by default. If you want to
#' customize the plot and only display specific variables, you can provide their
#' names as a character vector. For example, setting `varx = c("var1", "var2")`
#' will plot only "var1" and "var2". Note that the variable names should match
#' the column names in the dataset.
#' @param remove.zeros Specifies whether to remove zeros in the regression estimates
#' or shrinkage factors before calculating the averages and standard errors.
#' The default value is FALSE, indicating that zeros are not removed. If set to
#' TRUE, any zeros present in the regression estimates or shrinkage factors
#' will be excluded from the calculations. It is important to note that if
#' inclusion frequencies of variables are desired, zeros should not be removed.
#'@param plotype the plot to generate. Default is `plotype = "means"` which plots
#'the average values of a given metric. If distributions are required,
#'set `plotype = "distribution"` to generate boxplots.
#'@param legend.position the position of the legend in the plot. Default is to
#'place it on the right side of the plot
#'@param row,col One of `betatype` or `corrtype`,indicating the variables for
#'the rows and columns for plotting facet grid; note that row and col must be
#'different. Default is `row = betatype` and `col = corrtype`, so that the
#'plotting grid displays the given metric versus the SNR or sample size,
#'in a plotting grid with the coefficient types across the rows, and the
#'correlation levels across the columns.
#'@param linetype Specifies whether points with corresponding lines should be
#'plotted for the averages. The default is FALSE, which means only averages are
#'plotted. Please note that this option is only relevant when `plotype = "means"`.
#'@param add.se Specifies whether to include one standard error (SE) in
#'addition to the average values of the given metric. The default is FALSE. This
#'option is only relevant when `plotype = "means"`.
#' @param lwd,main,ylim,facet graphical parameters. see `[ggplot()]` and
#' `[plot]` for details
#' @inheritParams plot_metric
#' @export
plot_coefs <- function(x,
                       betatypes=c("a", "b", "c", "d"),
                       corrtypes=c("c1", "c2", "c3", "c4", "c5"),
                       samplesize = NULL,
                       snrx = NULL,
                       nvar = 15,
                       methods = NULL,
                       plotype = c("means", "distribution"),
                       legend.position=c("right", "bottom", "top", "left", "none"),
                       row=c("betatype", "corrtype"),
                       col=c("corrtype", "betatype"),
                       varx = NULL,
                       truncate.lower = NULL,
                       truncate.upper = NULL,
                       remove.zeros = FALSE,
                       linetype = FALSE,
                       add.se = FALSE,
                       lwd = 0.7,
                       points_shape = 16,
                       points_size = 2,
                       main = NULL,
                       ylim = NULL,
                       boxwidth = 0.5,
                       facet = FALSE){
  # Check for ggplot2 package
  if (!requireNamespace("ggplot2",quietly=TRUE))
    stop("Package ggplot2 not installed. Please install it!")

  # calculate the metric of interest if x is path where the rds files are stored
  #x <- combine_sim_coef(path, type = type)

  # assert that x must be of class "combinecoef"
  if (!inherits(x, "combinecoef"))
    stop("x must be an output from combine_sim_coef() function",call. = FALSE )

  # match arguments
  corrtypes <- match.arg(corrtypes)
  betatypes <- match.arg(betatypes)
  legend.position <- match.arg(legend.position)
  col<-match.arg(col)
  row<-match.arg(row)
  plotype <- match.arg(plotype)

  # Assert that row and col must be different
  if (row==col)
    stop("row and col must be different", call. = FALSE)

  # metric used in combine_sim_metrics()
  type <- attr(x, "metric")

  if(!corrtypes%in%unique(x$corrtype))
    stop("The corrtypes = ",corrtypes, " not used in the simulation", call. = FALSE)

  if(!betatypes%in%unique(x$betatype))
    stop("The betatypes = ",betatypes, " not used in the simulation", call. = FALSE)

  if(!is.null(samplesize) && !all(samplesize%in%unique(x$n)))
    stop("Not all sample sizes supplied were used in the simulation", call. = FALSE)

  if(!is.null(snrx) && !all(snrx%in%unique(x$snr)))
    stop("Not all SNR supplied were used in the simulation", call. = FALSE)

  # Set y-label and x-label
  ylab  = switch(type,
                 "beta" = "Regression Estimates",
                 "shrinkage"   = "Shrinkage Factors",
                 "inclusion" = "Inclusion Freq."
  )

  xlab <- "Variable"

  # all methods used in the analysis. If methods for plotting are
  #not provided we plot all methods
  methodx <- unique(x$method)

  if (!is.null(methods)) {
    invalid_methods <- methods[!methods %in% methodx]
    if (length(invalid_methods) > 0)
      stop("The methods ", invalid_methods, " are invalid", call. = FALSE)
    # use user methods
    methodx <- methods
  }

  # Deal with signal to noise ratio
  snr_used <- unique(x$snr)

  if (!is.null(snrx)) {
    snr_index <- which(!snrx%in%snr_used)
    if (length(snr_index)!=0)
      stop("The snr supplied was not used in the simulation. \n",
           sprintf("i This includes snr of: %s.",
                   paste0(snrx[snr_index], collapse = ", ")))
    # snr to be used for plotting
    snr_used <- unique(snrx)
  }

  # Deal with sample size
  n_used <- unique(x$n)
  if (!is.null(samplesize)) {
    ind <- which(!samplesize%in%n_used)
    if (length(ind)!=0)
      stop("The sample supplied was not used in the simulation. \n",
           sprintf("i This includes sample size of: %s.",
                   paste0(samplesize[ind], collapse = ", ")))
    n_used <- unique(samplesize)
  }

  # Deal number of covariates used in the model
  if (length(nvar)!=1)
    stop("! nvar must be a single number specifying the number of variables used.\n",
         sprintf("i The length of nvar is %d.", length(nvar)))

  nvar_used <- unique(x$nvars)
  nvar_index <- which(!nvar%in%nvar_used)

  if (length(nvar_index)!=0)
    stop("! nvar = ", nvar, " was not used in the simulation.", call. = FALSE)

  # Deal with the variables to be plotted
  var_used <- names(x)[-c(1:8)]

  if (!is.null(varx)) {
  index_invalid <- which(!varx%in%var_used)
    if (any(!varx%in%var_used))
      stop("Not all variables supplied in varx are subset of variables used in\n the simulation",
           sprintf("This applies to the following variables: %s.",
                   paste0(varx[index_invalid], collapse = ",")))
  # use the user input variables
  var_used <- varx
  }

  # Convert x to data.table to use its operation especially subsetting and
  # replacing 0s
  x <- data.table::as.data.table(x)
  x <- x[complete.cases(x)]

  # subset the data based on simulation parameters and convert from wide to long
  me1<- x[corrtype==corrtypes & betatype==betatypes &snr%in%snr_used & n%in%n_used & nvars==nvar & method%in%methodx]

  # Get rid of unnecessary variables as a result of different lengths
  me1 <- me1[,1:(8+nvar)]
  long_data <- data.table::melt(me1,
                                id.vars = c("scenario", "betatype","corrtype",
                                            "snr", "n","sigma", "nvars", "method"),
                                variable.name = "variable",
                                value.name = "metric")
  # subset long data based on variables: this is how is done in data.table
  long_data <- long_data[variable%in%var_used]
  #long_data <- data.table::as.data.table(dplyr::filter(long_data, variable%in%var_used))

  # Deal with truncation. Important for negative shrinkage factors
  if (type=="shrinkage"){
     if (!is.null(truncate.lower))
       long_data[metric <= truncate.lower, metric := truncate.lower]

    if (!is.null(truncate.upper))
      long_data[metric >= truncate.upper, metric := truncate.upper]
  }


   # Get rid of zeros in metric before aggregating. Just replace 0s by NA
  if (remove.zeros)long_data[metric == 0, metric := NA]
  # Remove rows with NA values
  long_data <- long_data[complete.cases(long_data)]
  # calculate the mean within a group with corresponding standard error and add
  # to the long_data by replicating the values within a group
  dfx <- long_data[, c(.SD,mean_value = mean(metric, na.rm = T),
                       se = std_mean(metric)), by = .(scenario,method, variable)]

  # rename betatype and corrtypes for plotting
  dfx[betatype == betatypes, betatype := paste0("Beta-type = ", betatypes)]
  dfx[corrtype == corrtypes, corrtype := paste0("Corr-type = ", corrtypes)]

  # all possible combination of snr and n when the length is different from 1
  sampx <- expand.grid(snr =snr_used, n=n_used)
  nx <- dim(sampx)[1L]

  # save plots as list
  plots <- vector("list", nx)

  if(plotype=="means"){
  # Get rid of the metric that has been replaced by the statistic.
  # Important for removing duplicates
    dfx <- subset(dfx, select = -metric)
    dfx <- dfx[!duplicated(dfx), ]

    for (v in 1:nx) {
        snx <-  sampx[v,1]
        nn <- sampx[v,2]
        datax <- dfx[snr == snx & n == nn]
        mainx <- paste0("n=", nn,":snr=",snx)
        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = factor(.data$variable),
                                                   y = .data$mean_value,
                                                   color = .data$method)) +
          ggplot2::geom_point(position = position_dodge(width = 0.5), shape = points_shape, size = points_size) +
          theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
        if (add.se) plots[[v]] <- plots[[v]] +
          geom_errorbar(aes(ymin = .data$mean_value - .data$se, ymax = .data$mean_value + .data$se),
                        width = 0.2, position = position_dodge(width = 0.5))
        if (linetype) plots[[v]] <- plots[[v]] + ggplot2::geom_line(aes(linetype = method, color = method, group = method),position = position_dodge(width = 0.5), lwd = lwd)
        if (!is.null(ylim)) plots[[v]] <- plots[[v]] + ylim(ylim)
        if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col)))
        if (is.null(main)) plots[[v]] <- plots[[v]] +
          ggplot2::ggtitle(paste0("snr=",snx, ":n = ", nn))

        }
  } else {
    # plot distributions
    for (v in 1:nx) {
      snx <-  sampx[v,1]
      nn <- sampx[v,2]
      datax <- dfx[snr == snx & n == nn]
      mainx <- paste0("n=", nn,":snr=",snx)
        plots[[v]] <- ggplot2::ggplot(data = datax,
                                      ggplot2::aes(x = .data$variable,
                                                   y = .data$metric,
                                                   colour = .data$method)) +
          ggplot2::geom_boxplot(width = boxwidth) +theme_bw() +  ggplot2::theme(legend.position = legend.position) +
          ggplot2::xlab(xlab) + ggplot2::ylab(ylab)

      if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col)))
      if (is.null(main)) plots[[v]] <- plots[[v]] + ggplot2::ggtitle(mainx)
      if (!is.null(ylim)) plots[[v]] <- plots[[v]] + ggplot2::coord_cartesian(ylim = ylim)
      if (facet) plots[[v]] <- plots[[v]] + ggplot2::facet_grid(formula(paste(row, "~", col)))
      if (is.null(main)) plots[[v]] <- plots[[v]] +
          ggplot2::ggtitle(paste0("snr=",snx, ":n = ", nn))

    }

  }
  plots
}

#' Plot regression estimate against their corresponding Shrinkage
#'
#' This function generates a plot to visualize the relationship between shrinkage
#' and regression estimates.
#'
#' @param beta1 Matrix of regression estimates from a variable selection approach.
#' @param beta2 Matrix of regression estimates from post-estimation shrinkage approach.
#' @param column Column that contains the variable of interest to be plotted.
#' @param remove.zeros Logical. Whether to remove rows with zero values in beta1 matrix.
#' @param points_shape Shape of the points in the plot (default: 16).
#' @param points_size Size of the points in the plot (default: 2).
#' @param betatrue True beta value for reference (default: 0).
#' @param linewidth Width of the dashed lines (default: 0.5).
#' @param truncate.lower A numeric value used to truncate negative shrinkage factors.
#' If a shrinkage factor is calculated to be less than this value, it will be set
#' to `truncate.lower`. Default is -1.
#' @param truncate.upper A numeric value used to truncate positive shrinkage factors.
#' If a shrinkage factor is calculated to be greater than this value, it will be set
#' to `truncate.upper`. Default is 1.25.
#' @return A plot visualizing the relationship between shrinkage and regression estimates.
#' @export
plot_coefs_shrinkage <- function(beta1, beta2, column = 1, remove.zeros = TRUE,
                                 points_shape = 16, points_size = 2, betatrue = 0,
                                 linewidth = 0.5, truncate.lower = -1, truncate.upper = 1.25) {

  # Check if beta1 and beta2 are matrices
  if (!is.matrix(beta1) || !is.matrix(beta2)) {
    stop("Both 'beta1' and 'beta2' must be matrices.")
  }

  # Check if column is a positive integer
  if (length(column) != 1 || !is.numeric(column) || column <= 0 || column %% 1 != 0) {
    stop("The 'column' must be a positive integer.", call. = FALSE)
  }

  # Check if dimensions of matrices are identical
  if (!identical(dim(beta1), dim(beta2))) {
    stop("Matrices have different dimensions.", call. = FALSE)
  }

  # Assert truncate is a numeric value
  if (!is.numeric(truncate.lower) || !is.numeric(truncate.upper)) {
    stop("Both truncate.lower and truncate.upper must be numeric values.", call. = FALSE)
  }


  # Subset the column of interest
  beta1 <- beta1[, column]
  beta2 <- beta2[, column]

  # Get rid of NA values in beta2
  index1 <- which(is.na(beta2))

  if (length(index1) != 0) {
    beta1 <- beta1[-index1]
    beta2 <- beta2[-index1]
  }

  # Optionally remove rows with zero values in beta1
  if (remove.zeros) {
    index2 <- which(beta1 != 0)
    if (length(index2) != 0) {
      shrinkage <- beta2[index2] / beta1[index2]
      dd <- data.frame(beta = beta1[index2], shrinkage = shrinkage)
    } else {
      shrinkage <- beta2 / beta1
      dd <- data.frame(beta = beta1, shrinkage = shrinkage)
    }
  } else {
    shrinkage <- beta2 / beta1
    shrinkage[is.nan(shrinkage)] <- 0
    dd <- data.frame(beta = beta1, shrinkage = shrinkage)
  }

  # Truncate
  dd$shrinkage <- ifelse(dd$shrinkage <= truncate.lower, truncate.lower, dd$shrinkage)
  dd$shrinkage <- ifelse(dd$shrinkage >= truncate.upper, truncate.upper, dd$shrinkage)
  #dd$shrinkage <- pmax(pmin(dd$shrinkage, truncate.upper), truncate.lower)


  # Create a signshrink column based on shrinkage values
  dd$signshrink <- ifelse(dd$shrinkage < 0, "neg", "pos")


  # Generate the plot
  ggplot(dd, aes(x = beta, y = shrinkage, color = signshrink))  +
    geom_point(shape = points_shape, size = points_size) +
    theme_bw() +
    ylab("Shrinkage") +
    xlab(expression(hat(beta)))  +
    geom_vline(xintercept = betatrue, linetype = "dashed", linewidth = linewidth) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = linewidth) +
    theme(legend.position = "none") + scale_color_manual(values = c("pos" = "black", "neg" = "red"))
}


#' Helper function used to calculate standard error of the mean
#'
#' @param x vector of numeric values. NA are removed before calculating the
#' statistic
#' @param na.rm a logical evaluating to TRUE or FALSE indicating whether NA values should be stripped
#' before computating standard deviation.
#' @return the standard error. if length of x is one we assume standard
#' deviation is 0
std_mean <- function(x, na.rm = TRUE){
  if (length(x)==1)
    return(0)
  sd(x, na.rm = na.rm) / sqrt(length(x))
}

#'@title Helper function that computes the sample median
#' @description
#' computes the sample median and a selected pair of outer quantiles
#' having equal tail areas.  It is identical to the function
#' \code{smedian.hilow()} in the \code{Hmisc} package.
#' @param x numeric vector from which NAs will be removed automatically
#' @param na.rm defaults to TRUE unlike built-in functions, so that by default
#' NAs are automatically removed
#' @param conf.int The coverage probability the outer quantiles should target.
#' When the default, 0.95, is used, the lower and upper quantiles computed are
#' 0.025 and 0.975.
#' @return returns the following summary statistics:
#' \item{Median}{Median value of x}
#' \item{Lower}{Lower quantile}
#' \item{Uppper}{Upper quantile}
smedian <- function (x, conf.int = 0.95, na.rm = TRUE) {
  quant <- quantile(x, probs = c(0.5, (1 - conf.int)/2, (1 + conf.int)/2),
                    na.rm = na.rm)
  names(quant) <- c("Median", "Lower", "Upper")
  quant
}

