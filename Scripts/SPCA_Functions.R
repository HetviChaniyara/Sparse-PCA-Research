# Cardinality and Equality Constrained PlS-SEM
# Hetvi Chaniyara
# Bachelor End Project
#Revised version December 2025 by Katrijn
#Focused on (checking) convergence and computational efficiency
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

#' CEC-PLS-SEM
#'
#' @param X data of size IxJ
#' @param R number of components
#' @param epsilon tolerance for convergence
#' @param phi number of nonzero weights
#' @param rho penalty tuning parameter
#' @param constrained 1/0 for constrained versus penalized setting (W=P)
#'
#' @returns Component weights and loadings
#' @export
#'
#' @examples
CEC_PLS_SEM <-function(X, R, epsilon, phi,rho, constrained, MaxIter){
  
  J = dim(X)[2] # number of columns
  I = dim(X)[1] # number of rows
  ssx <- sum(X^2)  #caching
  XtX <- t(X)%*%X  #caching
  iter <- 0
  convAO <- 0
  
  # Get initialized parameters
  params <- Initialize_parameters(X,R,phi)
  alpha <- params$alpha
  W <- params$W0
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
    message('Update P: Diff loss ', Lossc-LossuP)
    ####
    
    # Update weights
    eigenp <- eigen(t(P)%*%P)
    alpha_c <- alpha*eigenp$values[1] # learning step
    for (i in 1:4){#MM iterative procedure
      # Compute B
      B <- compute_B(X,W,P, alpha_c, XtX)
      # Compute W
      W <- compute_W_new(X, R, P, B, alpha_c, rho, U, phi)
      LossuW <- loss_function(X,W,P,rho,U)/ssx
      ####
      message('Update W: Diff loss ', LossuP-LossuW)
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
      message('Update U: Diff loss ', LossuW-LossuU)
      LossuP <- LossuU
      ####
    }
    
    #primary & secondary relative residuals
    r1 <- sum((W-P)^2)/sum(W^2)
    r2 <- sum((W-Wold)^2)/(sum(U^2)+1e-9)
    message('Primary relative residual:  ', r1)
    message('Secondary relative residual:  ', r2)
    
    # Calculate loss
    Lossu <- loss_function(X,W,P,rho,U)/ssx
    Lossvec <- c(Lossvec,Lossu)
    
    #Check for convergence or if maximum iterations are reached
    if (iter > MaxIter) {
      convAO <- 1
      cat("Maxiter")
    }
    
    # Relative Stopping Criterion
    relative_change <- (abs(Lossu - Lossc)) / abs(Lossc)
    
    if (relative_change < epsilon) {
      convAO <- 1
      cat("convergence")
    }
    
    print(paste("Iteration completed:", iter))
    iter <- iter + 1
    Lossc <- Lossu
  }
  
  results <- list('weights' = W, 'loadings' = P, 'Lossvec' = Lossvec, 'Residual' = Lossu, 'Scores'= T_scores, 'n_iterations'= iter)
  return(results)
}

########################################################################################################################################
# Helper Functions

Initialize_parameters <- function(X, R, phi) {
  
  J <- dim(X)[2] # number of columns
  I <- dim(X)[1] # number of rows
  svd_X <- svd(X)
  W_svd <- svd_X$v[, 1:R]
  alpha <- svd_X$d[1]^2 # max eigenvalue of X^TX, more efficient
  
  # Random components: note sum of sq. W from svd =1
  W_rand <- matrix(rnorm(length(W_svd), mean = 0, sd = 1/sqrt(J)), nrow = nrow(W_svd))
  
  # Weighted combination: 0.7 * SVD + 0.3 * random
  W0 <- 0.97*W_svd + 0.03*W_rand
  Wind <- order(abs(W0))#absolute value needed!
  W0[Wind[1:(J*R-phi)]] <- 0
  
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

compute_W_new <- function(X, R, P, B, alpha, rho, U, phi_prop) {
  
  W_new <- ((2 * alpha * B) + rho * (P - U)) / (2 * alpha + rho)
  K <- phi_prop
  n_total <- length(W_new)
  
  if (K < n_total) {
    # indices of entries sorted by increasing |W_new|
    impind <- order(abs(W_new), decreasing = FALSE)
    W_new[impind[1:(n_total - K)]] <- 0
  }
  
  return(W_new)
}

normalize_columns <- function(W, tol = 1e-10) {
  
  R <- ncol(W)
  
  for (r in 1:R) {
    norm_val <- sqrt(sum(W[, r]^2))
    
    if (norm_val < tol) {
      warning(paste("Column", r, 
                    "of W has (near) zero norm. Algorithm stopped."))
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
  
  return(list(precision = precision,recall = recall,f1 = f1_score,recovery = accuracy))
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



