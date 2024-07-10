#' Calculate the degrees of freedom for a model based on nonzero predictors.
#'
#' This function takes a matrix of predictors (\code{x}) for all variables and a 
#' vector of coefficients (\code{betas}), and calculates the degrees of freedom, 
#' which is the rank of the covariate matrix for the selected variables.
#'
#' @param x A matrix of predictors.
#' @param betas A matrix of coefficients for each value of lambda.
#'
#' @return The degrees of freedom for the model.
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 5)
#' x <- cbind(x, x[,2]+x[,4])
#' betas <- matrix(c(0, 0.5, 0, 1, 0, 1))
#' rownames(betas) = sprintf("x%d", 1:6)
#' calculate_df(x, betas)
#'
#' @export
calculate_df <- function(x, betas) {
  # Check if x is a matrix
  if (!is.matrix(x)) {
    stop("Input 'x' must be a matrix.", call. = FALSE)
  }
  
  # Check if betas is a numeric vector
  if (!is.matrix(betas)){
    stop("Input 'betas' must be a matrix.", call. = FALSE)
  }
  
  # Check if the number of columns in x matches the length of betas
  if (ncol(x) != nrow(betas)) {
    stop("Number of columns in 'x' must be equal to the column of 'betas'.", call. = FALSE)
  }
    # nonzero variables
    var_each_lambda <- create_list(betas)
    
    # set to NULL when no variable is selected 
    nvar_each_lambda <- lapply(var_each_lambda, function(v) length(v))
    empty_model <- which(nvar_each_lambda==0)
    
    if(length(empty_model)!=0){
      var_each_lambda[empty_model] <- list(NULL)
    }
    
    # calculate df for each lambda value using rank of the selected matrix of variables
    nv <- length(var_each_lambda)
    df <- numeric(nv)
    for(i in 1:nv){
    # select variables. null model is denoted by NULL
    x_selected <- x[, var_each_lambda[[i]], drop = FALSE]
    
    # Rank for selected variables
    df[i] <- qr(x_selected)$rank
    }

  return(df)
}

# create indices of variables selected
create_list <- function(mat){
  nv <- ncol(mat)
  nonzero <- vector(mode="list", length = nv)
  for(k in 1:nv){
    nonzero[[k]] <- which(mat[,k,drop = TRUE]!=0)
  }
  return(nonzero)
}

#' Compute the degrees of freedom for ridge regression.
#'
#' This function calculates the degrees of freedom for ridge regression for a given design matrix and a sequence of lambda values.
#'
#' @param x A numeric matrix representing the design matrix.
#' @param lambda A numeric vector representing the sequence of lambda values.
#'
#' @return A numeric vector containing the degrees of freedom for each lambda value.
#'
#' @examples
#' # Example usage:
#' x <- matrix(rnorm(20), ncol = 2)
#' lambda_values <- c(0.1, 0.01, 0.001)
#' df_ridge_result <- df_ridge(x, lambda_values)
#'
#' @references
#' see element of statistical learning for the formula. eqn 3.5
#'
#' @export
df_ridge <- function(x, lambda) {
  unlist(lapply(
    lambda,
    function(v) trace_matrix(x %*% solve(t(x) %*% x + v * diag(ncol(x))) %*% t(x))
  ))
}
#' Trace of a Matrix
#'
#' Calculates the trace of a given matrix.
#'
#' @param A The input matrix.
#' @return The sum of the diagonal elements of the matrix.
#' @details This function checks if the input is a matrix and returns the sum of its diagonal
#' elements.
#' @examples
#' A <- matrix(1:9, 3, 3)
#' trace_matrix(A)
#'
#' @export
trace_matrix <- function(A) { 
  if (!is.matrix(A)) {
    stop("A must be a matrix")
  }
  return(sum(diag(A)))
}


