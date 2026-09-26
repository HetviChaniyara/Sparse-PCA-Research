# Block GPower (l0) Simulation Benchmark
# September 2026
#
# Runs Scripts/GPower_Functions.R::gpower_block_cardinality() over the same
# simulated data grid as Demo/TSPCA_Sim.R (TSPCA) and Demo/SPCA_Sim.R
# (Zou's SPCA), with the same metric set, stability reporting and parallel
# execution -- see Demo/TSPCA_Sim.R's header for the full rationale behind both
# additions; only the parts specific to GPower are re-explained below.
#
# GPower has no multistart: gpower_block_cardinality() uses a deterministic
# SVD warm start (see Scripts/GPower_Functions.R), so unlike TSPCA
# there is exactly one fit per dataset, not N_MULTISTART restarts picking
# the best. That makes each job lighter than a TSPCA_Sim.R job, but the run is
# still worth parallelising: same per-dataset job granularity as TSPCA_Sim.R and
# SPCA_Sim.R, one job per (scheme, dataset-index) pair, flattened across
# all four schemes and distributed across workers up front.

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("../Scripts/TSPCA_Functions.R")
source("../Scripts/GPower_Functions.R")
library(dplyr)
library(doParallel)
library(foreach)

# ---- configuration ----------------------------------------------------------
# "total": cardinality budget shared across the whole JxR weight matrix
#   (matches the "block" convention TSPCA_Sim.R and seafar now run
#   under -- see Demo/TSPCA_Sim.R's header comment on CARDINALITY_TYPE).
# "per_component": each component gets its own fixed budget.
CARDINALITY_TYPE <- "total"
MU               <- 1     # equal per-component weighting, matching gPower's own default
MAX_ITER         <- 1000  # gpower_block_cardinality()'s own default
TOL              <- 1e-4  # gpower_block_cardinality()'s own default

N_CORES <- max(1, parallel::detectCores() - 1)

# ---- data-generation schemes (identical to Demo/TSPCA_Sim.R) --------------------
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

# ---- build a flat job list across ALL schemes -------------------------------
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
    source("../Scripts/GPower_Functions.R")
    library(dplyr)
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

      phi <- if (CARDINALITY_TYPE == "total") {
        round((1 - design$p_sparse[job$i]) * J * R)
      } else {
        round((1 - design$p_sparse[job$i]) * J) # recycled to every column by gpower_block_cardinality
      }

      res <- gpower_block_cardinality(X, R = R, phi = phi, mu = MU,
                                       cardinality_type = CARDINALITY_TYPE,
                                       max_iter = MAX_ITER, tol = TOL)

      W_aligned <- align_components(res$weights, trueW)
      P_aligned <- align_components(res$loadings, trueP)
      selection <- evaluate_variable_selection(trueW, W_aligned)
      bvm_W <- compute_bias_variance_mse(trueW, W_aligned)
      bvm_P <- compute_bias_variance_mse(trueP, P_aligned)
      bvm_W_P <- compute_bias_variance_mse(W_aligned, P_aligned)
      msd_W_P <- mean((W_aligned - P_aligned)^2)

      data.frame(
        Scheme = job$scheme,
        Folder = job$folder,
        Dataset = job$i,
        design[job$i, , drop = FALSE],
        Loss = NA, # gpower_block_cardinality() maximizes a different objective (sum of thresholded squared correlations), not a residual comparable to TSPCA's
        FEV = res$exp_var,
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
        Var_W = bvm_W$variance,
        Var_P = bvm_P$variance,
        msd_W_P = msd_W_P,
        W_Corr = diag(cor(W_aligned, trueW)) %>% mean(),
        P_Corr = diag(cor(P_aligned, trueP)) %>% mean(),
        Iterations = res$iterations,
        MSE_W_P = bvm_W_P$mse
      )
    }, error = function(e) {
      msg <- sprintf("[%s / %s%d] FAILED: %s", job$scheme, job$prefix, job$i, conditionMessage(e))
      cat(msg, "\n", file = "gpower_errors.log", append = TRUE)
      NULL
    })
  }

  # See Demo/TSPCA_Sim.R's matching comment: foreach's automatic export is
  # transitive (walks into try_job()'s body to find design_by_scheme, etc.),
  # so the explicit .export below is a redundant belt-and-suspenders
  # safeguard -- hence the suppressWarnings() around the harmless "already
  # exporting" warning.
  results_list <- suppressWarnings(
    foreach(job = jobs, .errorhandling = "pass",
            .export = c("try_job", "design_by_scheme", "CARDINALITY_TYPE",
                        "MU", "MAX_ITER", "TOL")) %dopar% {
      try_job(job)
    }
  )

  stopCluster(cl)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  cat(sprintf("Parallel run finished in %s mins across %d worker(s).\n", elapsed, N_CORES))

  is_bad <- vapply(results_list, function(x) is.null(x) || inherits(x, "error"), logical(1))
  n_failed <- sum(is_bad)
  if (n_failed > 0) {
    cat(sprintf("%d / %d jobs failed and were excluded (see gpower_errors.log for the ones caught inside try_job()).\n",
                 n_failed, length(jobs)))
  }
  results_list <- results_list[!is_bad]
} else {
  results_list <- list()
}

# ---- summarise and write results ------------------------------------------
cat("Compiling and writing results...\n")
if (length(results_list) == 0) {
  cat("No results were produced -- none of the configured data folders were found, or every job failed. Nothing written.\n")
} else {
  results_df <- dplyr::bind_rows(results_list)
  write.csv(results_df, "GPower_results.csv", row.names = FALSE)

  design_group_vars <- c("Scheme", "Folder", "n_variables", "s_size", "p_sparse", "n_components", "VAFx")
  design_group_vars <- intersect(design_group_vars, names(results_df))

  final_summary <- results_df %>%
    group_by(across(all_of(design_group_vars))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  write.csv(final_summary, "GPower_summary.csv", row.names = FALSE)

  stability_summary <- summarise_stability(results_df, group_vars = design_group_vars)
  write.csv(stability_summary, "GPower_stability.csv", row.names = FALSE)

  cat("Process finished successfully!\n")
}
