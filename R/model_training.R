#' Train Model with Knockoffs
#'
#' @param X A matrix of original features where columns represent features and rows represent samples
#' @param X_k A matrix of knockoff features
#' @param y Response variable
#' @param ... Additional parameters passed to the specific training method
#'
#' @return A list containing the trained model and other relevant information
#' @export
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @examples
#' # Train model with knockoffs
#' model <- model_training(X, X_k, y)

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
  return(as.numeric(abs(temp)))
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

# ---- Main Training Function ----

model_training <- function(X, X_k, y, ...) {
  # Create 5-fold cross-validation indices
  foldid <- sample(1:5, nrow(X), replace = TRUE)
  colnames(X) <- paste0("X", 1:ncol(X))
  colnames(X_k) <- paste0("X_k", 1:ncol(X_k))
  # Create main effects and squared terms
  fSX <- create.fSX(as.matrix(cbind(X, X_k)))
  # First-step Lasso
  cv.fit.1st <- glmnet::glmnet(fSX, y, standardize = TRUE, alpha = 1)
  beta_current <- cv.fit.1st$beta
  first_hit_50 <- 0
  tau <- NULL; order.index <- NULL
  for (iter_beta in 1:ncol(beta_current)) {
    if (sum(beta_current[, iter_beta] != 0) > 50) {
      beta <- as.matrix(beta_current[, iter_beta])
      T.all <- sapply(1:(ncol(X) * 2), get.T.j.fSX, x = cbind(X, X_k), beta = beta)
      T.all <- as.matrix(colMeans(T.all))
      T_0 <- as.matrix(abs(T.all[1:ncol(X)]))
      T_k <- matrix(abs(T.all[-(1:ncol(X))]), nrow = ncol(X), byrow = FALSE)
      MK.fit <- MK.statistic(T_0, T_k)
      tau <- MK.fit[, 'tau']
      order.index <- order(tau, decreasing = TRUE)
      if (sum(tau != 0) > 50) {
        first_hit_50 <- iter_beta
        break
      }
    }
  }
  # Fallback if no solution found
  if (first_hit_50 == 0) {
    first_hit_50 <- which.max(colSums(beta_current != 0))
  }
  # Refit Lasso with selected lambda
  cv.fit.first <- glmnet::glmnet(fSX, y, standardize = TRUE, alpha = 1, lambda = cv.fit.1st$lambda[first_hit_50])
  N <- sum(tau != 0)
  index.2nd <- intersect(order.index[1:N], which(tau != 0))
  index.2nd_k <- index.2nd
  # Create interaction terms for selected features
  fIX <- create.fIX(as.matrix(cbind(X[, index.2nd, drop = FALSE], X_k[, index.2nd_k, drop = FALSE])))
  # Second-step Lasso
  cv.fit.2nd <- glmnet::glmnet(cbind(fSX, fIX), y, standardize = TRUE, alpha = 1)
  prediction <- cbind(fSX, fIX) %*% cv.fit.2nd$beta
  y_replicate <- matrix(y, nrow = length(y), ncol = ncol(X), byrow = FALSE)
  BIC_vec <- (nrow(X) * log(colMeans((prediction - y_replicate)^2))) + ((colSums(cv.fit.2nd$beta != 0)) * log(nrow(X)))
  N_BIC <- which.min(BIC_vec)
  beta_current <- cv.fit.2nd$beta[-(1:ncol(fSX)), , drop = FALSE]
  beta.1st <- as.matrix(beta_current[, N_BIC])
  non_zero_rows <- which(beta.1st != 0)
  pairs_current <- rownames(beta.1st)[non_zero_rows]
  final_selected_set_k_biggest <- c()
  names_all <- colnames(fIX)
  final_interaction_added <- c()
  for (pair_i in seq_along(pairs_current)) {
    names_split <- strsplit(pairs_current[pair_i], ":")
    num_1 <- gsub("\\D", "", names_split[[1]][1])
    num_2 <- gsub("\\D", "", names_split[[1]][2])
    xi_xj <- paste0("X", num_1, ":X", num_2)
    xi_xkj <- paste0("X", num_1, ":X_k", num_2)
    xki_xj <- paste0("X", num_2, ":X_k", num_1)
    xki_xkj <- paste0("X_k", num_1, ":X_k", num_2)
    if (!(xi_xj %in% names_all)) xi_xj <- paste0("X", num_2, ":X", num_1)
    if (!(xi_xkj %in% names_all)) xi_xkj <- paste0("X_k", num_2, ":X", num_1)
    if (!(xki_xj %in% names_all)) xki_xj <- paste0("X_k", num_1, ":X", num_2)
    if (!(xki_xkj %in% names_all)) xki_xkj <- paste0("X_k", num_2, ":X_k", num_1)
    final_selected_set_k_biggest <- c(final_selected_set_k_biggest, xi_xj, xi_xkj, xki_xj, xki_xkj)
    final_interaction_added <- c(final_interaction_added, xi_xj)
  }
  final_interaction_added <- unique(final_interaction_added)
  final_selected_set_k_biggest <- unique(final_selected_set_k_biggest)
  fIX_filter <- fIX[, colnames(fIX) %in% final_selected_set_k_biggest]
  # Final cross-validated Lasso
  cv.fit.2nd_final <- glmnet::cv.glmnet(cbind(fSX, fIX_filter), y, standardize = TRUE, alpha = 1, foldid = foldid)
  lasso.fit <- glmnet::glmnet(cbind(fSX, fIX_filter), y, standardize = TRUE, alpha = 1, lambda = cv.fit.2nd_final$lambda[which.min(cv.fit.2nd_final$cvm)])
  return(list(
    trained_model = lasso.fit
  ))
} 
