# Tied Sparse PCA (TSPCA): sparse weights tied to the loadings (W = P)
# Hetvi Chaniyara
# Bachelor End Project
#Revised version December 2025 by Katrijn
#Focused on (checking) convergence and computational efficiency
#Revised further September 2026: cardinality can now be specified either as a
#single budget for the whole W matrix ("total", the original behaviour) or
#separately per component ("per_component"), and a matching-boundary bug in
#the initializer was fixed (see apply_cardinality()).
#Implementation allows two variants:
#1. Constrained approach
#min(W,P) ||X-XWP'||^2 s.t. Card(W) = K and P=W
#2. Penalized approach
#min(W,P) ||X-XWP'||^2 + rho/2||W-P||^2  s.t. Card(W) = K
#This is achieved by reformulating the objective to
#min(W,P,(U)) ||X-XWP'||^2 + rho/2||W-P+U||^2-rho/2||U||^2 s.t. Card(W) = K
#and fixing U to 0 or not

library(MASS)
library(gtools)

#############

#' TSPCA
#'
#' @param X data of size IxJ
#' @param R number of components
#' @param epsilon tolerance for convergence
#' @param phi number of nonzero weights. Either a single number, or (when
#'   cardinality_type = "per_component") a vector of length R giving a
#'   separate budget per component. When cardinality_type = "total", phi is
#'   the total number of nonzero entries allowed across the whole JxR weight
#'   matrix (the original behaviour); when cardinality_type = "per_component",
#'   a single phi is recycled to every column, or a length-R vector sets a
#'   different cardinality per component.
#' @param rho penalty tuning parameter
#' @param constrained 1/0 for constrained versus penalized setting (W=P)
#' @param verbose TRUE/FALSE for printing info iterates or not
#' @param cardinality_type "total" (default, matches the original behaviour)
#'   or "per_component". "total" ranks all J*R entries of W together and
#'   keeps the phi largest in magnitude, which can (and, for unbalanced data,
#'   sometimes does) allocate very unevenly across components -- including,
#'   in the worst case, a component with none, which normalize_columns()
#'   will stop on. "per_component" instead applies the cardinality budget
#'   within each column of W separately, guaranteeing every component keeps
#'   exactly phi (or phi[r]) nonzero weights.
#'
#' @returns Component weights and loadings
#' @export
#'
#' @examples
TSPCA <-function(X, INIT=NULL, R, epsilon, phi,rho, constrained, MaxIter,verbose=F,
                  cardinality_type = c("total", "per_component")){

  cardinality_type <- match.arg(cardinality_type)

  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows

  if (cardinality_type == "total" && length(phi) != 1) {
    stop("`phi` must be a single number when cardinality_type = 'total'.")
  }
  if (cardinality_type == "per_component" && !(length(phi) %in% c(1, R))) {
    stop("`phi` must have length 1 or R when cardinality_type = 'per_component'.")
  }

  ssx <- sum(X^2)  #caching
  XtX <- t(X)%*%X  #caching
  iter <- 0
  convAO <- 0

  # Get initialized parameters

  params <- Initialize_parameters(X,R,phi,cardinality_type)
  alpha <- params$alpha
  if (is.null(INIT)){
    W <- params$W0
  } else {
    W <- INIT
  }
  if (constrained==1){
    U <- params$U
  } else {
    U <- matrix(0,nrow = J, ncol = R)
  }
  
  # Initialize matrices and lists
  T_scores <- matrix(nrow = I, ncol = R)
  Lossc <- 1
  Lossvec <- Lossc
  
  # Update Loop
  while (convAO == 0) {
    Wold <- W #Wold needed for secondary residual
    
    # Update component scores
    T_scores <- X%*%W
    
    # Update loadings
    P = compute_P_new(X,W,T_scores,U,rho,R)
    LossuP <- loss_function(X,W,P,rho,U)/ssx
    ####
    if (verbose){
    message('Update P: Diff loss ', Lossc-LossuP)}
    ####
    
    # Update weights
    eigenp <- eigen(t(P)%*%P)
    alpha_c <- alpha*eigenp$values[1] # learning step
    for (i in 1:4){#MM iterative procedure
      # Compute B
      B <- compute_B(X,W,P, alpha_c, XtX)
      # Compute W
      W <- compute_W_new(X, R, P, B, alpha_c, rho, U, phi, cardinality_type)
      LossuW <- loss_function(X,W,P,rho,U)/ssx
      ####
      if (verbose){
        message('Update W: Diff loss ', LossuP-LossuW)}
      LossuP <- LossuW
      ####
    }
    
    norm_res <- normalize_columns(W)
    W <- norm_res$W
    
    if (norm_res$zero_column) {
      stop("Algorithm terminated due to zero column in W.")
    }
    
    # Update scaled variable
    if (constrained==1){
      U <- compute_U(U, W, P, rho)
      ####
      LossuU <- loss_function(X,W,P,rho,U)/ssx
      if (verbose){
      message('Update U: Diff loss ', LossuW-LossuU)}
      LossuP <- LossuU
      ####
    }
    
    #primary & secondary relative residuals
    r1 <- sum((W-P)^2)/sum(W^2)
    r2 <- sum((W-Wold)^2)/(sum(U^2)+1e-9)
    if (verbose){
    message('Primary relative residual:  ', r1)
    message('Secondary relative residual:  ', r2)
    }
    # Calculate loss
    Lossu <- loss_function(X,W,P,rho,U)/ssx
    Lossvec <- c(Lossvec,Lossu)
    
    #Check for convergence or if maximum iterations are reached
    if (iter > MaxIter) {
      convAO <- 1
      message("Maxiter")
    }
    
    # Relative Stopping Criterion
    relative_change <- (abs(Lossu - Lossc)) / abs(Lossc)
    
    if (relative_change < epsilon) {
      convAO <- 1
      message("convergence")
    }
    
    if (verbose){
    print(paste("Iteration completed:", iter))}
    iter <- iter + 1
    Lossc <- Lossu
  }
  
  results <- list('weights' = W, 'loadings' = P, 'Lossvec' = Lossvec, 'Residual' = Lossu, 'Scores'= T_scores, 'n_iterations'= iter)
  return(results)
}

########################################################################################################################################
# Helper Functions

#' Apply a cardinality (hard-thresholding) constraint to a weight matrix
#'
#' Zeroes out all but the largest-magnitude entries of M, either counting
#' across the whole matrix ("total") or separately within each column
#' ("per_component"). Used by both Initialize_parameters() and
#' compute_W_new() so the two stay consistent, and guards the boundary case
#' (phi >= number of available entries) that a bare `1:(n-phi)` used to get
#' wrong in R when n-phi <= 0.
#'
#' @param M a JxR matrix of candidate weights (or scores to threshold)
#' @param phi total cardinality (scalar) when cardinality_type = "total";
#'   scalar (recycled to every column) or length-R vector of per-column
#'   cardinalities when cardinality_type = "per_component"
#' @param cardinality_type "total" or "per_component"
#' @returns M with all but the retained entries set to 0
apply_cardinality <- function(M, phi, cardinality_type = c("total", "per_component")) {

  cardinality_type <- match.arg(cardinality_type)
  J <- nrow(M)
  R <- ncol(M)

  if (cardinality_type == "total") {

    if (length(phi) != 1) stop("`phi` must be a single number when cardinality_type = 'total'.")
    K <- phi
    n_total <- J * R
    if (K >= n_total) return(M) # nothing to drop; guards the old off-by-one edge case
    drop_idx <- order(abs(M), decreasing = FALSE)[seq_len(n_total - K)]
    M[drop_idx] <- 0
    return(M)

  }

  # per_component: enforce the budget within each column separately
  if (!(length(phi) %in% c(1, R))) {
    stop("`phi` must have length 1 or R when cardinality_type = 'per_component'.")
  }
  phi_r <- if (length(phi) == 1) rep(phi, R) else phi

  for (r in seq_len(R)) {
    Kr <- phi_r[r]
    if (Kr < J) { # guards the same off-by-one edge case per column
      drop_idx <- order(abs(M[, r]), decreasing = FALSE)[seq_len(J - Kr)]
      M[drop_idx, r] <- 0
    }
  }
  return(M)
}

Initialize_parameters <- function(X, R, phi, cardinality_type = c("total", "per_component")) {

  cardinality_type <- match.arg(cardinality_type)
  J <- dim(X)[2] # number of columns
  I <- dim(X)[1] # number of rows
  svd_X <- svd(X,R,R)
  #random part
  W_rand <- matrix(rnorm(J * R, sd=1/sqrt(J)), ncol = R, nrow = J)
  #rational part
  W_rat_unr <- svd_X$v # %*% diag(svd_X$d[1:R]) / sqrt(I)
  if (R > 1) {  #rotation to simple structure
    varimax_res <- stats::varimax(W_rat_unr, normalize = FALSE)
    W_rat <- W_rat_unr %*% varimax_res$rotmat
  }
  #W_rat <- W_rat_unr
  #W_svd <- svd_X$v[, 1:R]
  alpha <- svd_X$d[1]^2 # max eigenvalue of X^TX, more efficient

  # Random components: note sum of sq. W from svd =1
  #W_rand <- matrix(rnorm(length(W_svd), mean = 0, sd = 1/sqrt(J)), nrow = nrow(W_svd))

  # Weighted combination: 0.7 * SVD + 0.3 * random
  W0 <- 0.8*W_rat + 0.2*W_rand
  W0 <- apply_cardinality(W0, phi, cardinality_type)

  U <- matrix(0, nrow = J, ncol = R) # Initialize to 0

  return(list(W0 = W0, U = U, alpha = alpha))
}

compute_P_new <- function(X, W, T_scores, U, rho, R) {
  
  # Calculate X^T XW
  XtXW <- t(X) %*% T_scores
  
  # Add regularization term rho * (W + U)
  regularization_term <- rho * (W + U)
  
  # Combine the terms
  term1 <- 2 * XtXW + regularization_term
  
  # Calculate (2 * W^T X^T X W + rho * I)
  I <- diag(R) 
  term2 <- 2 *(t(T_scores) %*% T_scores) + (rho * I)#! factor 2
  
  # Inverse of term2
  term2_inv <- solve(term2)#instead of ginv as this is a well defined problem
  
  # Multiply term1 by the inverse of term2
  P_new <- term1 %*% term2_inv
  
  return(P_new)
}

compute_B <- function(X,W,P, alpha,XTX){
  # Compute: PX_kron^T*PX_kron*vec(W) by identity = vec(X^TXWP_TP)
  term1 = (XTX %*% W %*% t(P) %*% P)
  
  # PX_kron^T *vec(X)
  term2 = (XTX %*% P)
  
  # Subtract term2 from term 1 and dividing by alpha
  term3 = term1 - term2
  term4 = term3/alpha
  
  # Subtract vec_W - term 4
  B = W - term4
  
  return(B)
}

# compute_W_new <- function(X, R, P, B, alpha, rho, U, phi_prop) {
# 
#  W_new <- ((2 * alpha * B) + rho * (P - U)) / (2 * alpha + rho)
#  # Coefficients with smallest bjr^2 + (Ujr-Pjr)^2 set to 0
#  term1 <- alpha*(B^2)
#  term2 <- 0.5*rho*((U-P)^2)
#  impind <- order(term1+term2,decreasing = FALSE)
#  J <- dim(X)[2]
#  W_new[impind[1:(J*R-phi_prop)]] <- 0
# 
#  return(W_new)
# }

compute_W_new <- function(X, R, P, B, alpha, rho, U, phi_prop, cardinality_type = c("total", "per_component")) {

  cardinality_type <- match.arg(cardinality_type)
  W_new <- ((2 * alpha * B) + rho * (P - U)) / (2 * alpha + rho)
  W_new <- apply_cardinality(W_new, phi_prop, cardinality_type)

  return(W_new)
}

normalize_columns <- function(W, tol = 1e-10) {
  
  R <- ncol(W)
  
  for (r in 1:R) {
    norm_val <- sqrt(sum(W[, r]^2))
    
    if (norm_val < tol) {
      warning(paste("Column", r, 
                    "of W has zero norm. Algorithm stopped."))
      return(list(W = W, zero_column = TRUE))
    }
    
    W[, r] <- W[, r] / norm_val
  }
  
  return(list(W = W, zero_column = FALSE))
}

compute_U <- function(U,W,P,rho){
  
  # Update U - without rho
  U_new <- U + (W- P)
  
  return(U_new)
}

loss_function <-function(X,W,P,rho,U){
  # Loss function
  term1 <- sum((X - X %*% W %*% t(P))^2)
  term2 <- (rho/2)*sum((W-P+U)^2)
  term3 <- (rho/2)*sum(U^2)
  total_loss <- term1+term2-term3
  return(total_loss)
}

###############################################################################################################################
# Evaluation Metrics Functions

#' Evaluate estimated vs. true zero/nonzero (variable-selection) status
#'
#' Revised September 2026: now also returns the raw confusion-matrix counts
#' and the false positive / false negative RATES (FPR, FNR) individually,
#' not just the blended precision/recall/f1/accuracy. In a simulation study
#' where every replication redraws its own random true support, "stability
#' of zero/nonzero status" can't be assessed positionally (item j has no
#' fixed identity across replications, unlike the Big5 real-data illustration
#' where the same 240 items are fixed across resamples) -- but the
#' *replication-to-replication variability of these rates*, computed here
#' per replication and then summarised (mean + SD/IQR) across the
#' n_replications datasets sharing a design cell in the calling script, is a
#' meaningful and directly interpretable stand-in: SD(FPR) across
#' replications reflects how consistently the method avoids spuriously
#' selecting truly-zero entries ("stability of the zero status"); SD(FNR)
#' reflects how consistently it retains truly-nonzero entries ("stability of
#' the nonzero status"). All previously-returned fields are unchanged, so
#' this is a backward-compatible extension.
#'
#' @param W_true JxR true weight/loading matrix (or a sub-block of one)
#' @param W_estimated JxR estimated weight/loading matrix, same shape,
#'   already aligned to W_true's column order (e.g. via align_components())
#' @returns list with precision, recall, f1, recovery (accuracy) as before,
#'   plus TP, FP, FN, TN, fpr (= FP / (FP+TN)), fnr (= FN / (FN+TP))
evaluate_variable_selection <- function(W_true, W_estimated) {

  # Checking which and how many coefficients are exactly 0
  W_true_bin <- ifelse(W_true != 0, 1, 0)
  W_est_bin <- ifelse(W_estimated != 0, 1, 0)

  TP <- sum(W_true_bin == 1 & W_est_bin == 1)
  FP <- sum(W_true_bin == 0 & W_est_bin == 1)
  FN <- sum(W_true_bin == 1 & W_est_bin == 0)
  TN <- sum(W_true_bin == 0 & W_est_bin == 0)

  precision <- TP / (TP + FP + 1e-8)
  recall <- TP / (TP + FN + 1e-8)
  f1_score <- 2 * (precision * recall) / (precision + recall + 1e-8)
  accuracy <- (TP + TN) / (TP + FP + FN + TN)
  fpr <- FP / (FP + TN + 1e-8)  # rate of truly-zero entries wrongly kept nonzero
  fnr <- FN / (FN + TP + 1e-8)  # rate of truly-nonzero entries wrongly zeroed

  return(list(precision = precision, recall = recall, f1 = f1_score, recovery = accuracy,
              TP = TP, FP = FP, FN = FN, TN = TN, fpr = fpr, fnr = fnr))
}

#' Summarise between-replication stability of zero/nonzero selection metrics
#'
#' For a simulation grid where each replication within a design cell draws
#' its OWN random true support (so positional zero/nonzero-flip comparison
#' across replications, as done for the fixed real-data Big5 illustration,
#' isn't meaningful -- see evaluate_variable_selection()'s docs above), this
#' computes, per design cell, the mean AND the between-replication dispersion
#' (SD and IQR) of a set of selection-quality metrics. Low dispersion across
#' the independently-generated replications of a cell means the method
#' reliably achieves similar (accurate or inaccurate) zero/nonzero
#' classification every time; high dispersion means performance is volatile
#' from dataset to dataset even though the design (sample size, sparsity,
#' etc.) is held fixed.
#'
#' @param results_df a data.frame with one row per replication (as produced
#'   by the main simulation loop), containing `group_vars` plus numeric
#'   columns named in `metric_vars`
#' @param group_vars character vector of column names identifying a design
#'   cell (e.g. c("Scheme","n_variables","s_size","p_sparse","n_components","VAFx"))
#' @param metric_vars character vector of numeric column names to summarise;
#'   defaults to the selection-stability-relevant set
#' @returns a data.frame, one row per design cell, with `<metric>_mean`,
#'   `<metric>_sd`, and `<metric>_iqr` for every entry in metric_vars, plus
#'   `n_replications` (how many replications contributed, after dropping NA)
summarise_stability <- function(results_df, group_vars,
                                 metric_vars = c("Recovery_Rate", "Precision", "Recall",
                                                  "F1", "FPR", "FNR")) {
  metric_vars <- intersect(metric_vars, names(results_df))
  results_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      n_replications = sum(!is.na(.data[[metric_vars[1]]])),
      dplyr::across(
        dplyr::all_of(metric_vars),
        list(mean = ~mean(.x, na.rm = TRUE),
             sd   = ~stats::sd(.x, na.rm = TRUE),
             iqr  = ~stats::IQR(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

reconstruction_metrics <- function(X, W, P) {
  
  # Reconstruction Metrics
  X_hat <- X %*% W %*% t(P)
  error_matrix <- X - X_hat
  mse <- mean(error_matrix^2)
  var_explained <- 1 - (sum(error_matrix^2) / sum((X - mean(X))^2))
  
  return(list(mse = mse, R2 = var_explained))
}

score_metrics <- function(est, true) {
  
  # General Function For MAE, RMSE and Corrleation
  mae <- mean(abs(est - true))
  rmse <- sqrt(mean((est - true)^2))
  corrs <- diag(cor(est, true)) # assumes same column order
  avg_corr <- mean(corrs)
  
  return(list(mae = mae, rmse = rmse, correlation = avg_corr))
}

align_components <- function(est, true) {
  # Try combinations to see which estimated composite is corresponding one in the true matrix
  n_comp <- ncol(true)
  perm <- permutations(n_comp, n_comp)
  
  best_perm <- NULL
  best_score <- Inf
  
  # Selects permutation with best score and returns that order of composites
  for (i in 1:nrow(perm)) {
    aligned_est <- est[, perm[i, ]]
    
    for (j in 1:n_comp) {
      correlation_val <- cor(aligned_est[, j], true[, j])
      
      # check if corr is NA 
      if (!is.na(correlation_val) && correlation_val < 0) {
        aligned_est[, j] <- -aligned_est[, j]
      }
    }
    
    score <- sum((aligned_est - true)^2)
    if (score < best_score) {
      best_score <- score
      best_perm <- aligned_est
    }
  }
  return(best_perm)
}

#' Align a weight/loading matrix to a reference zero/nonzero pattern
#'
#' Sparse PCA / TSPCA solutions carry no inherent component ordering across
#' independent runs (e.g. resamples, or different random starts): the
#' column that recovers "trait 1" in one run may end up in position 3 in
#' another. Comparing zero/nonzero status column-by-column across such runs
#' without first re-ordering columns will (mis)count real instability from
#' apparent instability caused by permutation alone. This finds, by brute
#' force over all column permutations (only feasible for small n_comp, as
#' used here with R = 5), the permutation whose nonzero pattern agrees most
#' with `ref_pattern`, and returns `est` reordered accordingly.
#'
#' @param est a JxR matrix of estimated weights or loadings
#' @param ref_pattern a JxR indicator matrix (1 = nonzero, 0 = zero) giving
#'   the expected zero/nonzero pattern that `est`'s columns should be
#'   matched against
#'
#' @returns list(perm = permutation applied to est's columns, aligned = est
#'   with columns reordered by perm)
align_to_zero_pattern <- function(est, ref_pattern) {
  n_comp <- ncol(ref_pattern)
  est_pattern <- (est != 0) * 1
  perm <- permutations(n_comp, n_comp)

  best_perm <- perm[1, ]
  best_score <- -Inf
  for (i in 1:nrow(perm)) {
    agreement <- sum(est_pattern[, perm[i, ]] == ref_pattern)
    if (agreement > best_score) {
      best_score <- agreement
      best_perm <- perm[i, ]
    }
  }
  return(list(perm = best_perm, aligned = est[, best_perm]))
}

compute_bias_variance_mse <- function(W_true, W_est) {
  
  # Compute bias, variance and MSE
  W_true_vec <- as.vector(W_true)
  W_est_vec <- as.vector(W_est)
  bias <- mean(W_est_vec - W_true_vec)
  variance <- var(W_est_vec - W_true_vec)
  mse <- mean((W_est_vec - W_true_vec)^2)
  
  return(list(bias = bias, variance = variance, mse = mse))
}

sparsity_level <- function(W) {
  
  # Checks the sparsity of the parameter
  total_elements <- length(W)
  zero_elements <- sum(W == 0)
  return(zero_elements / total_elements)
}

compute_vaf <- function(X, W, P) {
  # Variance Accounted For calculation
  X_hat <- X %*% W %*% t(P)
  sum_sq_error <- sum((X - X_hat)^2)
  total_variance <- sum(X^2)
  vaf <- 1 - (sum_sq_error / total_variance)
  return(vaf)
}



