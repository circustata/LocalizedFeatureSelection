#' Perform Localized Feature Selection
#'
#' This function implements the full localized feature selection pipeline:
#' 1. Use a trained lasso model (with interactions)
#' 2. Compute the augmented selection matrix and individual q-values for evaluation data
#'
#' @param trained_model Trained lasso model (from model_training)
#' @param X A matrix of original features (training or evaluation data)
#' @param X_k A matrix of knockoff features (corresponding to X)
#' @param fdr_threshold Target FDR threshold for selection (default 0.1)
#'
#' @return A list containing:
#'   \itemize{
#'     \item augmented_selection_matrix: The normalized selection matrix as defined in the paper
#'     \item Q: q-value matrix, Q[i, j] <= FDR targeted threshold means the j-th feature is selected for the i-th individual
#'     \item W_out: Localized W statistics (individualized feature importance)
#'   }
#' @export
#'
#' @examples
#' result <- feature_selection(trained_model, X, X_k, fdr_threshold = 0.1)

# ---- Helper Functions ----

# Create main effects and squared terms
create.fSX <- function(x) {
  x <- data.frame(x)
  temp.x2 <- x^2
  colnames(temp.x2) <- paste0(colnames(x), 'squared')
  temp <- cbind(x, temp.x2)
  return(as.matrix(temp))
}

# Create interaction terms (excluding self-interactions)
create.fIX <- function(x) {
  x <- data.frame(x)
  temp <- model.matrix(~ .^2, data = x)[, -1, drop = FALSE]
  final_columns <- c()
  column_names <- colnames(temp)
  for (col_name in column_names) {
    if (grepl(":", col_name)) {
      parts <- strsplit(col_name, ":")[[1]]
      numbers <- gsub("\\D", "", parts)
      if (length(unique(numbers)) > 1) {
        final_columns <- c(final_columns, col_name)
      }
    }
  }
  return(as.matrix(temp[, final_columns, drop = FALSE]))
}

# Calculate feature statistics for main effects
get.T.j.fSX <- function(j, x, beta) {
  x <- t(x)
  current_name <- rownames(x)[j]
  temp <- beta[which(rownames(beta) == current_name)]
  temp <- temp + 2 * beta[which(rownames(beta) == paste0(current_name, 'squared'))] * x[j, , drop = FALSE]
  return(as.numeric(temp))
}

# Calculate feature statistics for interactions
get.T.j.fIX <- function(j, x, beta) {
  x <- t(x)
  current_name <- rownames(x)[j]
  temp <- 0 * x[j, , drop = FALSE]
  regex <- paste0("(^|:)", current_name, "($|:)")
  indices <- grep(regex, rownames(beta), value = TRUE)
  if (length(indices) > 0) {
    for (names in indices) {
      other_part <- gsub(paste0(current_name, ":|:", current_name), "", names)
      temp <- temp + beta[names, ] * x[other_part, , drop = FALSE]
    }
  }
  return(as.numeric(temp))
}

# Compute knockoff statistics
MK.statistic <- function(T_0, T_k, method = 'median') {
  T_0 <- as.matrix(T_0)
  T_k <- as.matrix(T_k)
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0
  which.max.alt <- function(x) {
    temp.index <- which(x == max(x))
    if (length(temp.index) != 1) {
      return(temp.index[2])
    } else {
      return(temp.index[1])
    }
  }
  kappa <- apply(T.temp, 1, which.max.alt) - 1
  if (method == 'max') {
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, sort(x, decreasing = TRUE)[2])
  }
  if (method == 'median') {
    Get.OtherMedian <- function(x) median(x[-which.max(x)])
    tau <- apply(T.temp, 1, max) - apply(T.temp, 1, Get.OtherMedian)
  }
  return(cbind(kappa, tau))
}

# Compute q-values for knockoff statistics
MK.q.byStat <- function(kappa, tau, M = 1, Rej.Bound = 10000) {
  b <- order(tau, decreasing = TRUE)
  c_0 <- kappa[b] == 0
  ratio <- c(); temp_0 <- 0
  for (i in 1:length(b)) {
    temp_0 <- temp_0 + c_0[i]
    temp_1 <- i - temp_0
    temp_ratio <- (1/M + 1/M * temp_1) / max(1, temp_0)
    ratio <- c(ratio, temp_ratio)
    if (i > Rej.Bound) { break }
  }
  q <- rep(1, length(tau))
  if (length(which(tau[b] > 0)) != 0) {
    index_bound <- max(which(tau[b] > 0))
    for (i in 1:length(b)) {
      temp.index <- i:min(length(b), Rej.Bound, index_bound)
      if (length(temp.index) == 0) { next }
      q[b[i]] <- min(ratio[temp.index]) * c_0[i] + 1 - c_0[i]
      if (i > Rej.Bound) { break }
    }
    q[q > 1] <- 1
  }
  return(q)
}

# ---- Main Feature Selection Function ----

feature_selection <- function(trained_model, X, X_k, fdr_threshold = 0.1) {
  # Extract model coefficients
  beta <- trained_model$beta
  fSX <- create.fSX(as.matrix(cbind(X, X_k)))
  # Main effects and interactions
  beta.1st <- beta[(1:ncol(fSX)), , drop = FALSE]
  T.1st <- sapply(1:(ncol(X) * 2), get.T.j.fSX, x = cbind(X, X_k), beta = beta.1st)
  beta.2nd <- beta[-(1:ncol(fSX)), , drop = FALSE]
  fIX <- create.fIX(as.matrix(cbind(X, X_k)))
  T.2nd <- sapply(1:(ncol(X) * 2), get.T.j.fIX, x = cbind(X, X_k), beta = beta.2nd)
  T_all.abs <- abs(T.1st + T.2nd)
  T_0.all <- T_all.abs[, 1:ncol(X), drop = FALSE]
  T_k.all <- T_all.abs[, -(1:ncol(X)), drop = FALSE]

  n_eval <- nrow(X)
  Q <- matrix(NA, nrow = n_eval, ncol = ncol(X))
  W_out <- matrix(NA, nrow = n_eval, ncol = ncol(X))
  for (k in 1:n_eval) {
    T_0 <- as.matrix(T_0.all[k, ])
    T_k <- matrix(T_k.all[k, ], nrow = ncol(X), byrow = FALSE)
    MK.stat <- MK.statistic(T_0, T_k)
    W_out[k, ] <- ((MK.stat[, 'kappa'] == 0) - (MK.stat[, 'kappa'] == 1)) * MK.stat[, 'tau']
    q <- MK.q.byStat(MK.stat[, 'kappa'], MK.stat[, 'tau'], ncol(T_k), Rej.Bound = 10000)
    Q[k, ] <- q
  }

  # Compute augmented selection matrix: normalized W_out for Q < threshold, else 0
  augmented_selection_matrix <- matrix(0, nrow = n_eval, ncol = ncol(X))
  for (j in 1:ncol(X)) {
    max_wj <- max(abs(W_out[, j]))
    if (max_wj == 0) max_wj <- 1
    for (i in 1:n_eval) {
      if (Q[i, j] <= fdr_threshold) {
        augmented_selection_matrix[i, j] <- W_out[i, j] / max_wj
      }
    }
  }
  colnames(augmented_selection_matrix) <- colnames(X)

  return(list(
    augmented_selection_matrix = augmented_selection_matrix,
    Q = Q,
    W_out = W_out
  ))
} 

