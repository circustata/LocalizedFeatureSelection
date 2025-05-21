#' Generate Knockoffs using IPAD procedure
#'
#' @param X A matrix of original features (rows: samples, columns: features)
#' @param threshold Threshold for determining hyperparameter r (default: 0.01)
#'
#' @return A list containing:
#'   \itemize{
#'     \item X_k: Matrix of knockoff features
#'     \item X: Original feature matrix
#'   }
#' @export
#'
#' @examples
#' # Generate knockoffs for the data
#' result <- knockoffs_generate(X, threshold = 0.01)
knockoffs_generate <- function(X, threshold = 0.01) {
  n <- nrow(X)
  p <- ncol(X)

  # Eigen decomposition of X %*% t(X)
  a <- X %*% t(X)
  eigen_fit <- eigen(a)

  # Select r: the smallest r such that the average off-diagonal of residual covariance < threshold
  for (r in 1:p) {
    F_tilde <- sqrt(n) * eigen_fit$vectors[, 1:r]
    C <- F_tilde %*% (t(F_tilde) %*% X) / n
    E <- X - C
    Var_E <- cov(E)
    off_diag <- (sum(abs(Var_E)) - sum(diag(abs(Var_E)))) / (p * (p - 1))
    if (off_diag < threshold) {
      break
    }
  }

  # Generate knockoff features
  F_tilde <- sqrt(n) * eigen_fit$vectors[, 1:r]
  C <- F_tilde %*% (t(F_tilde) %*% X) / n
  E <- X - C
  E_e <- apply(E, 2, sd)
  E_k <- t(E_e * t(matrix(rnorm(n * p, 0, 1), n, p)))
  X_k <- C + E_k

  return(list(
    X_k = X_k,
    X = X
  ))
} 