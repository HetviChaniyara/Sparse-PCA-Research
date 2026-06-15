# Zou's Sparse PCA (Elastic Net) Benchmark - Dynamic Prefix Mode
# Updated March 2026

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("../Scripts/SPCA_Functions.R")
library(dplyr)
library(elasticnet) 

# Add all folders you want to process
folders <- c("../Scripts/DATA-R-W-Sparse")

for (f in folders) {
  # Load the design parameters for the data
  # Checking both common variations of the typo in "simulation"
  info_file <- file.path(f, "Info_simulation.RData")
  if (!file.exists(info_file)) info_file <- file.path(f, "Info_simulaiton.RData")
  
  if (!file.exists(info_file)) {
    cat("Skipping folder - missing simulation info file in:", f, "\n")
    next
  }
  
  load(info_file)
  # Handle case sensitivity variations in the object name safely
  design <- if(exists("Info_simulation")) Info_simulation$design_matrix_replication else Infor_simulation$design_matrix_replication
  total_files <- nrow(design) 
  
  # DYNAMIC PREFIX DETECTION: Determine the prefix based on folder name
  if (grepl("WP-Sparse", f, ignore.case = TRUE)) {
    prefix <- "WPsparse"
  } else if (grepl("P-Sparse", f, ignore.case = TRUE)) {
    prefix <- "Psparse"
  } else {
    prefix <- "Wsparse" # Default fallback to W-Sparse
  }
  
  benchmark_list <- vector("list", total_files)
  cat("Starting benchmark for folder:", f, "using prefix:", prefix, "for", total_files, "files...\n")
  
  start_time <- Sys.time()
  
  for (i in 1:total_files) {
    data_file <- file.path(f, paste0(prefix, i, ".RData"))
    if (!file.exists(data_file)) {
      next
    }
    
    tryCatch({
      load(data_file)
      
      X <- out$X
      R <- out$k
      
      # DYNAMIC MATRIX EXTRACTION
      # If it's a WP dataset, the true sparse weights are saved as 'W2'
      if (prefix == "WPsparse") {
        true_W <- out$W2  
        true_P <- out$P   
      } else {
        # Standard W-sparse or P-sparse files use 'W' and 'P'
        true_W <- out$W  
        true_P <- out$P  
      }
      
      rm(out) 
      
      phi_val <- colSums(true_W != 0) 
      
      # Fit Elastic Net
      enet_fit <- elasticnet::spca(X, K = R, para = phi_val, sparse = "varnum")
      
      # Component Alignment
      W_aligned <- align_components(enet_fit$loadings, true_W)
      rm(enet_fit) 
      
      # Fast Regression for P_aligned
      Z_scores  <- X %*% W_aligned
      tZ_Z      <- crossprod(Z_scores)
      tX_Z      <- crossprod(X, Z_scores)
      P_aligned <- tX_Z %*% solve(tZ_Z)
      
      # Compute Analytics Metrics
      selection <- evaluate_variable_selection(true_W, W_aligned)
      bvm_W     <- compute_bias_variance_mse(true_W, W_aligned)
      bvm_P     <- compute_bias_variance_mse(true_P, P_aligned) 
      bvm_W_P   <- compute_bias_variance_mse(W_aligned, P_aligned)
      msd_W_P <- mean((W_aligned - P_aligned)^2)
      
      w_corrs <- diag(cor(W_aligned, true_W))
      w_corrs[is.na(w_corrs)] <- 0
      p_corrs <- diag(cor(P_aligned, true_P))
      p_corrs[is.na(p_corrs)] <- 0
      
      # Save results row
      benchmark_list[[i]] <- data.frame(
        Method = "Zou_SPCA_ENet",
        Folder = f, 
        Dataset = i, 
        design[i, , drop = FALSE], # Keeps the exact 'n_components' column structure
        Loss = NA, 
        FEV = compute_vaf(X, W_aligned, P_aligned),
        Recovery_Rate = selection$recovery,
        MSE_W = bvm_W$mse,
        MSE_P = bvm_P$mse,
        Bias_W = bvm_W$bias,
        Bias_P = bvm_P$bias,
        Var_W = bvm_W$variance,
        Var_P = bvm_P$variance,
        W_Corr = mean(w_corrs),
        P_Corr = mean(p_corrs), 
        Iterations = NA, 
        MSE_W_P = bvm_W_P$mse,
        msd_W_P = msd_W_P
      )
      
    }, error = function(e) {
      cat("\nError reading dataset index", i, ":", conditionMessage(e), "\n")
    })
    
    rm(X, true_W, true_P, W_aligned, P_aligned, Z_scores)
    
    if (i %% 200 == 0) {
      elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
      cat(sprintf("Folder: %s | Progress: %d / %d | Elapsed: %s mins\n", 
                  prefix, i, total_files, elapsed))
      gc(verbose = FALSE) 
    }
  }
  
  cat("Compiling and summarizing datasets...\n")
  
  benchmark_summary <- bind_rows(benchmark_list) %>%
    group_by(Folder, n_variables, s_size, p_sparse, n_components, VAFx) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
  
  output_name <- paste0("ElasticNet_Summary_", prefix, ".csv")
  write.csv(benchmark_summary, output_name, row.names = FALSE)
  cat("Benchmarking finished successfully! Saved to:", output_name, "\n")
}