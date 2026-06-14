# Zou's Sparse PCA (Elastic Net) Benchmark
# Updated March 2026 - Hetvi Chaniyara
# Matches CEC-PLS-SEM loop and metric structure

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("../Scripts/SPCA_Functions.R")
library(dplyr)
library(elasticnet) 

# list of all data folders 
folders <- c("../Scripts/DATA-R-W-Sparse", "../Scripts/DATA-R-P-Sparse")
benchmark_list <- list()

for (f in folders) {
  # load the design parameters for the data
  load(file.path(f, "Info_simulation.RData"))
  design <- Info_simulation$design_matrix_replication
  prefix <- ifelse(grepl("W-Sparse", f), "Wsparse", "Psparse")
  
  for (i in 1:180) {
    data_file <- file.path(f, paste0(prefix, i, ".RData"))
    if(!file.exists(data_file)) {
      cat("File missing:", data_file, "\n")
      next
    }
    load(data_file)
    
    X <- out$X
    R <- out$k
    J <- ncol(X)
    
    phi_val <- colSums(out$W != 0) 
    
    # uses SVD to start so no multistart needed
    enet_fit <- elasticnet::spca(X, K = R, para = phi_val, type = "predictor", sparse = "varnum")
    
    # aligned components
    W_aligned <- align_components(enet_fit$loadings, out$W)
    
    # Calculate optimal P_aligned for the benchmark
    # P = (X'Z)(Z'Z)^-1
    ####
    #Z_scores <- X %*% W_aligned
    #P_aligned <- t(X) %*% Z_scores %*% solve(t(Z_scores) %*% Z_scores)
    ##### ADAPTED TO ORIGINAL ELASTICNET P
    alpha <- X%*%t(X)%*%W_aligned
    z <- svd(alpha)
    P_aligned <- (z$u) %*% t(z$v)
    
    # metrics
    selection <- evaluate_variable_selection(out$W, W_aligned)
    bvm_W <- compute_bias_variance_mse(out$W, W_aligned)
    bvm_P <- compute_bias_variance_mse(out$P, P_aligned) 
    bvm_W_P <- compute_bias_variance_mse(W_aligned, P_aligned)
    w_corrs <- diag(cor(W_aligned, out$W))
    w_corrs[is.na(w_corrs)] <- 0
    p_corrs <- diag(cor(P_aligned, out$P))
    p_corrs[is.na(p_corrs)] <- 0
    
    benchmark_list[[length(benchmark_list)+1]] <- data.frame(
      Method = "Zou_SPCA_ENet",
      Folder = f, 
      Dataset = i, 
      design[i,],
      Loss = NA, # not same residual like CEC-PLS-SEM
      VAF = compute_vaf(X, W_aligned, P_aligned),
      Recovery_Rate = selection$recovery,
      MSE_W = bvm_W$mse,
      MSE_P = bvm_P$mse,
      Bias_W = bvm_W$bias,
      Bias_P = bvm_P$bias,
      W_Corr = mean(w_corrs),
      P_Corr = mean(p_corrs), 
      Iterations = NA, # no iterations recorded from method
      MSE_W_P = bvm_W_P$mse
    )
    cat("Folder:", f, "Dataset:", i, "Benchmark Complete\n")
  }
}

# summarise results
benchmark_summary <- do.call(rbind, benchmark_list) %>%
  group_by(Folder, n_variables, s_size, p_sparse, VAFx) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

write.csv(benchmark_summary, "all_elastic_net.csv", row.names = FALSE)