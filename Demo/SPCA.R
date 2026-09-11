# TSPCA Research
# Updated Version March 2026 Hetvi Chaniyara
# Runs TSPCA for various data folders
# Incorporates the changes of Katrijn's December 2025 version
#
# Revised September 2026: now runs automatically over all data-generation
# schemes (W-sparse, P-sparse, WP-sparse) without hand-editing which `out$`
# field to read -- see the `conditions` table below, which encodes exactly
# what used to be manual "# Change to W2 if WPsparse" / "X <- out$X2"
# swapping. Also fixes: the write.csv() call used to receive the raw list
# of per-dataset data.frames instead of the row-bound table (write.csv() on
# a list does not do what you want -- it silently mis-serializes); the
# `if (prefix==1)` branch below was always FALSE (prefix is a string, never
# the number 1) so variable-selection recovery was always evaluated against
# the wrong ground truth for every condition; and folders/files that don't
# exist yet (e.g. simulations still being generated) are now skipped with a
# clear message instead of crashing the whole run.

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

source("../Scripts/SPCA_Functions.R")
library(dplyr)

# ---- configuration --------------------------------------------------------
# "total": cardinality budget is shared across the whole JxR weight matrix
#   (the original behaviour -- can allocate unevenly across components).
# "per_component": each component gets its own budget, derived the same way
#   Elastic_net.R already derives its per-column oracle target. Set this to
#   "per_component" for a fairer, crash-safer comparison; see
#   SPCA_Functions.R::apply_cardinality for the mechanics.
CARDINALITY_TYPE <- "total"
CONSTRAINED      <- 0     # 0 = penalized (W minus P), 1 = equality-constrained (W = P)
N_MULTISTART     <- 10
MAX_ITER         <- 100

# ---- data-generation schemes ----------------------------------------------
# Each entry says: which folder/file prefix to read, and how to pull out the
# analysis data (X) and ground truth (W, P) from the saved `out` list for
# that scheme. The WP-sparse folder yields *two* conditions from the same
# files: the asymmetric branch (X / W / P, loadings sparse, weights free)
# and the symmetric branch (X2 / W2, the literal W = P case) -- see the note
# at the top of Scripts/WPsparse_Data.R.
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

# The three generator scripts don't agree on the info-file name or the name
# of the object saved inside it ("Info_simulaiton.RData" / Infor_simulation
# for the W- and WP-sparse generators, vs "Info_simulation.RData" /
# Info_simulation for the P-sparse one) -- handle both so callers don't need
# to know which.
load_design <- function(folder) {
  candidates <- c("Info_simulaiton.RData", "Info_simulation.RData")
  found <- candidates[file.exists(file.path(folder, candidates))]
  if (length(found) == 0) {
    stop("No Info_simulation(.RData) file found in ", folder)
  }
  env <- new.env()
  load(file.path(folder, found[1]), envir = env)
  obj_name <- intersect(c("Infor_simulation", "Info_simulation"), ls(env))
  if (length(obj_name) == 0) {
    stop("Loaded ", found[1], " from ", folder,
         " but found neither Infor_simulation nor Info_simulation inside.")
  }
  get(obj_name[1], envir = env)$design_matrix_replication
}

results_list <- list()

for (cond in conditions) {

  f <- cond$folder

  if (!dir.exists(f)) {
    cat(sprintf("Skipping scheme '%s': folder %s does not exist (not generated on this machine yet).\n",
                cond$scheme, f))
    next
  }

  design <- load_design(f)
  n_datasets <- nrow(design)

  start_time <- Sys.time()
  n_missing <- 0

  for (i in seq_len(n_datasets)) {

    data_file <- file.path(f, paste0(cond$prefix, i, ".RData"))
    if (!file.exists(data_file)) {
      n_missing <- n_missing + 1
      next
    }
    load(data_file) # provides `out`

    true_params <- cond$extract(out)
    X <- true_params$X
    trueW <- true_params$W
    trueP <- true_params$P

    R <- out$k
    J <- ncol(X)

    phi <- if (CARDINALITY_TYPE == "total") {
      round((1 - design$p_sparse[i]) * J * R)
    } else {
      round((1 - design$p_sparse[i]) * J) # recycled to every column by CEC_PLS_SEM
    }
    rho <- sum(X^2) / R

    # run with multistart
    best_res <- NULL; best_loss <- Inf
    for (m in 1:N_MULTISTART) {
      set.seed(100 + m)
      res <- CEC_PLS_SEM(X, INIT = NULL, R, 1e-8, phi, rho, constrained = CONSTRAINED,
                          MaxIter = MAX_ITER, cardinality_type = CARDINALITY_TYPE)
      if (res$Residual < best_loss) { best_loss <- res$Residual; best_res <- res }
    }

    # metrics calculation -- always against this condition's own ground truth,
    # no more "Change to W2 if WPsparse" hand-editing needed
    W_aligned <- align_components(best_res$weights, trueW)
    P_aligned <- align_components(best_res$loadings, trueP)
    selection <- evaluate_variable_selection(trueW, W_aligned)
    bvm_W <- compute_bias_variance_mse(trueW, W_aligned)
    bvm_P <- compute_bias_variance_mse(trueP, P_aligned)
    bvm_W_P <- compute_bias_variance_mse(W_aligned, P_aligned)
    msd_W_P <- mean((W_aligned - P_aligned)^2)

    # Storing results in the dataframe
    results_list[[length(results_list) + 1]] <- data.frame(
      Scheme = cond$scheme,
      Folder = f,
      Dataset = i,
      design[i, , drop = FALSE],
      Loss = best_res$Residual,
      FEV = compute_vaf(X, best_res$weights, best_res$loadings),
      Recovery_Rate = selection$recovery,
      MSE_W = bvm_W$mse,
      MSE_P = bvm_P$mse,
      Bias_W = bvm_W$bias,
      Bias_P = bvm_P$bias,
      Var_W = bvm_W$variance,
      Var_P = bvm_P$variance,
      msd_W_P = msd_W_P,
      W_Corr = diag(cor(W_aligned, trueW)) %>% mean(),
      P_Corr = diag(cor(P_aligned, trueP)) %>% mean(),
      Iterations = best_res$n_iterations,
      MSE_W_P = bvm_W_P$mse
    )

    # flush memory
    rm(X, trueW, trueP, true_params, W_aligned, P_aligned, best_res, selection,
       bvm_W, bvm_P, bvm_W_P, msd_W_P, out, res)

    # garbage collection / progress
    if (i %% 100 == 0) {
      elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
      cat(sprintf("Scheme: %s | Completed: %d / %d | Elapsed: %s mins\n",
                   cond$scheme, i, n_datasets, elapsed))
      gc(verbose = FALSE)
    }
  }

  if (n_missing > 0) {
    cat(sprintf("Scheme '%s': %d / %d data files were missing and skipped (generation likely incomplete).\n",
                 cond$scheme, n_missing, n_datasets))
  }
}

# summarise and write results
cat("Compiling and writing results...\n")
if (length(results_list) == 0) {
  cat("No results were produced -- none of the configured data folders were found. Nothing written.\n")
} else {
  results_df <- do.call(rbind, results_list)
  write.csv(results_df, "CEC_PLS_SEM_results.csv", row.names = FALSE)

  final_summary <- results_df %>%
    group_by(Scheme, Folder, n_variables, s_size, p_sparse, VAFx) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  write.csv(final_summary, "CEC_PLS_SEM_summary.csv", row.names = FALSE)

  cat("Process finished successfully!\n")
}
