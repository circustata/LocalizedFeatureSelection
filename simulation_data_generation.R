#' Generate Simulation Data (Scenario 3)
#'
#' This function generates simulation data according to Scenario 3 in our paper.
#' It returns both the main data (X_original, Z1, Z2, y) and the evaluation data (out_X_original, out_Z1, out_Z2), as well as the combined X and out_X.
#' @param n_train Number of training samples (default 1000)
#' @param n_eval Number of evaluation samples (default 2000)
#' @param p Number of variables/features (default 100)
#' @param k Number of nonzero coefficients (default 20)
#' @param rho Correlation parameter for Toeplitz covariance (default 0.5)
#' @return A list with:
#'   \item{X_original}{Main feature matrix for training (n_train x p)}
#'   \item{Z1}{Covariate matrix 1 for training (n_train x 2)}
#'   \item{Z2}{Covariate matrix 2 for training (n_train x 1)}
#'   \item{X}{Combined training data matrix (n_train x (p+3)), with features and covariates}
#'   \item{y}{Response vector for training data}
#'   \item{out_X_original}{Main feature matrix for evaluation (n_eval x p)}
#'   \item{out_Z1}{Covariate matrix 1 for evaluation (n_eval x 2)}
#'   \item{out_Z2}{Covariate matrix 2 for evaluation (n_eval x 1)}
#'   \item{out_X}{Combined evaluation data matrix (n_eval x (p+3)), with features and covariates}
#'   \item{beta_true}{True regression coefficients for the first p features}
#'   \item{nonzero}{Indices (1:p) of true signal features among the first p columns}
#'   \item{covariate_index}{Indices of covariate columns (Z1, Z2) in X and out_X}
#' @export
generate_simulation_data <- function(
  n_train = 1000,
  n_eval = 2000,
  p = 100,
  k = 20,
  rho = 0.5
) {
  set.seed(123) # For reproducibility, you can remove or set as argument

  # Covariance structure
  mu <- rep(0, p)
  Sigma <- toeplitz(rho^(0:(p-1)))
  Sigma_all <- matrix(0, p+3, p+3)
  Sigma_all[1:p, 1:p] <- Sigma
  diag(Sigma_all) <- 1

  # Nonzero coefficients (true signals among first p features)
  nonzero <- sort(sample(1:p, k))
  amplitude <- 12 / sqrt(n_train)
  beta_true <- amplitude * (1:p %in% nonzero)

  # Covariate columns (Z1, Z2) are at the end
  covariate_index <- (p + 1):(p + 3)

  # Training data
  X_original <- matrix(rnorm(n_train * p), n_train) %*% chol(Sigma)
  Z1 <- matrix(rbinom(n_train * 2, size = 1, prob = 0.5), n_train)
  Z2 <- matrix(runif(n_train, min = 0.2, max = 1), n_train)
  X <- cbind(X_original, Z1, Z2)
  y <- as.matrix(
    X_original %*% beta_true + rnorm(n_train)
  )

  # Evaluation data
  out_X_original <- matrix(rnorm(n_eval * p), n_eval) %*% chol(Sigma)
  out_Z1 <- matrix(rbinom(n_eval * 2, size = 1, prob = 0.5), n_eval)
  out_Z2 <- matrix(runif(n_eval, min = 0.2, max = 1), n_eval)
  out_X <- cbind(out_X_original, out_Z1, out_Z2)

  return(list(
    X_original = X_original,
    Z1 = Z1,
    Z2 = Z2,
    X = X,
    y = y,
    out_X_original = out_X_original,
    out_Z1 = out_Z1,
    out_Z2 = out_Z2,
    out_X = out_X,
    beta_true = beta_true,
    nonzero = nonzero,
    covariate_index = covariate_index
  ))
} 