# Zou's Sparse PCA (Elastic Net) Benchmark
# Updated March 2026 - Hetvi Chaniyara
# Matches TSPCA_Sim.R loop and metric structure
#
# Revised September 2026: same two additions made to Demo/TSPCA_Sim.R, applied
# here so the benchmark stays comparable.
#  1. Stability of zero/nonzero status. See Demo/TSPCA_Sim.R's header for the
#     full rationale -- in short, each replication draws its own random
#     true support, so between-replication SD/IQR of the selection-quality
#     metrics (not positional zero/nonzero agreement) is what's meaningful.
#     Written to Elastic_Net_stability.csv, one row per design cell.
#  2. Parallel execution under doParallel/foreach, one job per dataset
#     (elasticnet::spca has no multistart to also parallelise, unlike
#     TSPCA, so the per-job unit of work is lighter here -- but still
#     worth distributing across cores over hundreds of datasets).
#  Also: folder existence / missing-file handling now matches TSPCA_Sim.R
#  (skips cleanly with a message instead of failing on the whole run), and
#  the info file is read with the same load_design() as TSPCA_Sim.R.
#
# Revised September 2026 (second pass): added the WP-sparse folder, split
# into its two conditions the same way Demo/TSPCA_Sim.R and Demo/GPower_Sim.R
# already do -- the asymmetric branch (X/W/P, loadings sparse, weights
# free) and the symmetric branch (X2/W2, the literal W = P case; see the
# note at the top of Scripts/WPsparse_Data.R). Conditions are now specified
# via an `extract` function per scheme (matching the other three scripts)
# instead of reading out$X/out$W/out$P directly, so this loop drives all
# four schemes uniformly.

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("../Scripts/TSPCA_Functions.R")
library(dplyr)
library(elasticnet)
library(doParallel)
library(foreach)

# Number of parallel worker processes -- same convention as TSPCA_Sim.R.
N_CORES <- max(1, parallel::detectCores() - 1)

# ---- data-generation schemes (identical to Demo/TSPCA_Sim.R / Demo/GPower_Sim.R) -
conditions <- list(
  list(folder = "../Scripts/DATA-R-W-Sparse",  prefix = "Wsparse",  scheme = "W_sparse",
       extract = function(out) list(X = out$X,  W = out$W,  P = out$P)),
  list(folder = "../Scripts/DATA-R-P-Sparse",  prefix = "Psparse",  scheme = "P_sparse",
       extract = function(out) list(X = out$X,  W = out$W,  P = out$P)),
  list(folder = "../Scripts/DATA-R-WP-Sparse", prefix = "WPsparse", scheme = "WP_sparse_asymmetric",
       extract = function(out) list(X = out$X,  W = out$W,  P = out$P)),
  list(folder = "../Scripts/DATA-R-WP-Sparse", prefix = "WPsparse", scheme = "WP_sparse_symmetric_WeqP",
       extract = function(out) list(X = out$X2, W = out$W2, P = out$W2))
)

# Reads the design matrix that the data-generation scripts save alongside the
# datasets in Info_simulation.RData.
load_design <- function(folder) {
  env <- new.env()
  load(file.path(folder, "Info_simulation.RData"), envir = env)
  env$Info_simulation$design_matrix_replication
}

# ---- build a flat job list across both schemes ------------------------------
jobs <- list()
design_by_scheme <- list()

for (cond in conditions) {
  f <- cond$folder
  if (!dir.exists(f)) {
    cat(sprintf("Skipping scheme '%s': folder %s does not exist (not generated on this machine yet).\n",
                cond$scheme, f))
    next
  }
  design <- load_design(f)
  design_by_scheme[[cond$scheme]] <- design
  n_datasets <- nrow(design)
  n_missing <- 0
  for (i in seq_len(n_datasets)) {
    data_file <- file.path(f, paste0(cond$prefix, i, ".RData"))
    if (!file.exists(data_file)) {
      n_missing <- n_missing + 1
      next
    }
    jobs[[length(jobs) + 1]] <- list(scheme = cond$scheme, folder = f, prefix = cond$prefix,
                                      i = i, data_file = data_file, extract = cond$extract)
  }
  if (n_missing > 0) {
    cat(sprintf("Scheme '%s': %d / %d data files were missing and will be skipped (generation likely incomplete).\n",
                 cond$scheme, n_missing, n_datasets))
  }
}

cat(sprintf("Total jobs to run: %d, across %d worker(s).\n", length(jobs), N_CORES))

# ---- parallel execution -----------------------------------------------------
if (length(jobs) > 0) {

  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    source("../Scripts/TSPCA_Functions.R")
    library(dplyr)
    library(elasticnet)
  })

  start_time <- Sys.time()

  try_job <- function(job) {
    tryCatch({
      design <- design_by_scheme[[job$scheme]]
      load(job$data_file) # provides `out`

      true_params <- job$extract(out)
      X <- true_params$X
      trueW <- true_params$W
      trueP <- true_params$P
      R <- out$k
      J <- ncol(X)

      # elasticnet::spca is given the TRUE per-component cardinality as an
      # oracle target (it needs a per-column budget, unlike TSPCA's
      # shared/total option) -- unchanged from the original benchmark.
      phi_val <- colSums(trueW != 0)

      # SVD-based, deterministic -- no multistart needed (unchanged).
      enet_fit <- elasticnet::spca(X, K = R, para = phi_val, type = "predictor", sparse = "varnum")

      W_aligned <- align_components(enet_fit$loadings, trueW)

      # Optimal P given W, via the original elasticnet-style polar solution
      # (unchanged from the prior version).
      alpha <- t(X) %*% X %*% W_aligned
      z <- svd(alpha)
      P_aligned <- (z$u) %*% t(z$v)

      selection <- evaluate_variable_selection(trueW, W_aligned)
      bvm_W <- compute_bias_variance_mse(trueW, W_aligned)
      bvm_P <- compute_bias_variance_mse(trueP, P_aligned)
      bvm_W_P <- compute_bias_variance_mse(W_aligned, P_aligned)
      w_corrs <- diag(cor(W_aligned, trueW))
      w_corrs[is.na(w_corrs)] <- 0
      p_corrs <- diag(cor(P_aligned, trueP))
      p_corrs[is.na(p_corrs)] <- 0

      data.frame(
        Method = "Zou_SPCA_ENet",
        Scheme = job$scheme,
        Folder = job$folder,
        Dataset = job$i,
        design[job$i, , drop = FALSE],
        Loss = NA, # not the same residual as TSPCA
        VAF = compute_vaf(X, W_aligned, P_aligned),
        Recovery_Rate = selection$recovery,
        Precision = selection$precision,
        Recall = selection$recall,
        F1 = selection$f1,
        FPR = selection$fpr,
        FNR = selection$fnr,
        MSE_W = bvm_W$mse,
        MSE_P = bvm_P$mse,
        Bias_W = bvm_W$bias,
        Bias_P = bvm_P$bias,
        W_Corr = mean(w_corrs),
        P_Corr = mean(p_corrs),
        Iterations = NA, # no iterations recorded from this method
        MSE_W_P = bvm_W_P$mse
      )
    }, error = function(e) {
      msg <- sprintf("[%s / %s%d] FAILED: %s", job$scheme, job$prefix, job$i, conditionMessage(e))
      cat(msg, "\n", file = "elastic_net_errors.log", append = TRUE)
      NULL
    })
  }

  # See TSPCA_Sim.R's matching comment: foreach's automatic export is transitive
  # (walks into try_job()'s body to find design_by_scheme), so the explicit
  # .export below is a redundant belt-and-suspenders safeguard -- hence the
  # suppressWarnings() around the harmless "already exporting" warning.
  results_list <- suppressWarnings(
    foreach(job = jobs, .errorhandling = "pass",
            .export = c("try_job", "design_by_scheme")) %dopar% {
      try_job(job)
    }
  )

  stopCluster(cl)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  cat(sprintf("Parallel run finished in %s mins across %d worker(s).\n", elapsed, N_CORES))

  is_bad <- vapply(results_list, function(x) is.null(x) || inherits(x, "error"), logical(1))
  n_failed <- sum(is_bad)
  if (n_failed > 0) {
    cat(sprintf("%d / %d jobs failed and were excluded (see elastic_net_errors.log for the ones caught inside try_job()).\n",
                 n_failed, length(jobs)))
  }
  results_list <- results_list[!is_bad]
} else {
  results_list <- list()
}

# ---- summarise and write results --------------------------------------------
cat("Compiling and writing results...\n")
if (length(results_list) == 0) {
  cat("No results were produced -- none of the configured data folders were found, or every job failed. Nothing written.\n")
} else {
  results_df <- dplyr::bind_rows(results_list)
  write.csv(results_df, "Elastic_Net_results.csv", row.names = FALSE)

  design_group_vars <- c("Scheme", "Folder", "n_variables", "s_size", "p_sparse", "n_components", "VAFx")
  design_group_vars <- intersect(design_group_vars, names(results_df))

  benchmark_summary <- results_df %>%
    group_by(across(all_of(design_group_vars))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  write.csv(benchmark_summary, "all_elastic_net.csv", row.names = FALSE)

  # Stability of zero/nonzero selection status across replications within
  # each design cell -- see TSPCA_Sim.R's header and TSPCA_Functions.R::
  # summarise_stability() for what this does and does not capture.
  stability_summary <- summarise_stability(results_df, group_vars = design_group_vars)
  write.csv(stability_summary, "Elastic_Net_stability.csv", row.names = FALSE)

  cat("Process finished successfully!\n")
}
