# seafar (USLPCA) Simulation Benchmark
# September 2026
#
# Runs seafar::seafar_multistart() over the same simulated data grid as
# Demo/TSPCA_Sim.R (TSPCA), Demo/SPCA_Sim.R (Zou's SPCA) and
# Demo/GPower_Sim.R (block GPower), with the same metric set, stability
# reporting and parallel execution -- see Demo/TSPCA_Sim.R's header for the full
# rationale behind both additions; only what's specific to seafar is
# re-explained below.
#
# seafar's own multistart is internal (that's what the function's name
# means): seafar_multistart() already restarts internally and returns its
# best solution, so -- like GPower and unlike TSPCA -- there is no
# outer N_MULTISTART loop here.
#
# seafar has no separate weights/loadings distinction: seafar_multistart()
# returns a single sparse loading matrix ($loadings) that IS what scores
# are computed from ($scores = X %*% loadings, used directly -- see
# Illustration_Big5.qmd), unlike TSPCA/GPower/spca, which each return
# a weight matrix (used to compute scores) and a separate, generally denser
# reconstruction loading matrix. So below, the one matrix seafar returns is
# aligned once and evaluated against BOTH trueW and trueP separately (as
# the WP_sparse_symmetric_WeqP condition in TSPCA_Sim.R already does for the
# W = P case) -- for the W_sparse/P_sparse/WP_sparse_asymmetric conditions,
# where trueW != trueP, this means neither comparison should be expected to
# look as clean as it would for a method that fits the two separately.
#
# Function signature confirmed from Katrijn's own working usage in
# Demo/Illustration_Big5.qmd / .Rhistory:
#   seafar:::seafar_multistart(X, nfactors = R, C = <total nonzero budget>,
#                               INIT = "rational")
# C is a TOTAL cardinality across the whole loading matrix (e.g. 240 = 5
# components x 48 items/trait in the Big5 illustration) -- there is no
# documented per-component variant, so this script always computes phi the
# "total" way regardless of CARDINALITY_TYPE below (kept only so the
# column matches the other three scripts' output for easy comparison).
#
# NOTE: this script could not be executed against a live copy of the seafar
# package in the environment it was written in (GitHub access for
# remotes::install_github() was blocked there) -- it was written directly
# from the confirmed call signature above and smoke-tested with a stub
# standing in for seafar:::seafar_multistart() to verify the surrounding
# pipeline (job building, parallel execution, metric computation, stability
# aggregation) only. Please do a first small run (e.g. restrict `jobs` to a
# handful of datasets) before committing to the full grid, to confirm the
# real function's behavior/timing on your machine.

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

source("../Scripts/TSPCA_Functions.R")
library(dplyr)
library(seafar)
library(doParallel)
library(foreach)

# ---- configuration ----------------------------------------------------------
CARDINALITY_TYPE <- "total"  # seafar's C is always a total/shared budget; see note above
SEAFAR_INIT      <- "rational"

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
    library(dplyr)
    library(seafar)
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

      phi <- round((1 - design$p_sparse[job$i]) * J * R)  # total budget, see header note

      res <- seafar:::seafar_multistart(X, nfactors = R, C = phi, INIT = SEAFAR_INIT)

      # seafar returns one sparse matrix -- see header note. Align it once
      # and evaluate against both trueW and trueP.
      Loadings_aligned <- align_components(res$loadings, trueW)

      selection_W <- evaluate_variable_selection(trueW, Loadings_aligned)
      selection_P <- evaluate_variable_selection(trueP, Loadings_aligned)
      bvm_W <- compute_bias_variance_mse(trueW, Loadings_aligned)
      bvm_P <- compute_bias_variance_mse(trueP, Loadings_aligned)

      data.frame(
        Scheme = job$scheme,
        Folder = job$folder,
        Dataset = job$i,
        design[job$i, , drop = FALSE],
        Loss = NA, # seafar_multistart() has no residual directly comparable to TSPCA's
        FEV = compute_vaf(X, res$loadings, res$loadings),
        Recovery_Rate = selection_W$recovery,
        Precision = selection_W$precision,
        Recall = selection_W$recall,
        F1 = selection_W$f1,
        FPR = selection_W$fpr,
        FNR = selection_W$fnr,
        Recovery_Rate_vsP = selection_P$recovery,
        Precision_vsP = selection_P$precision,
        Recall_vsP = selection_P$recall,
        F1_vsP = selection_P$f1,
        FPR_vsP = selection_P$fpr,
        FNR_vsP = selection_P$fnr,
        MSE_W = bvm_W$mse,
        MSE_P = bvm_P$mse,
        Bias_W = bvm_W$bias,
        Bias_P = bvm_P$bias,
        Var_W = bvm_W$variance,
        Var_P = bvm_P$variance,
        W_Corr = diag(cor(Loadings_aligned, trueW)) %>% mean(),
        P_Corr = diag(cor(Loadings_aligned, trueP)) %>% mean(),
        Iterations = NA
      )
    }, error = function(e) {
      msg <- sprintf("[%s / %s%d] FAILED: %s", job$scheme, job$prefix, job$i, conditionMessage(e))
      cat(msg, "\n", file = "seafar_errors.log", append = TRUE)
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
            .export = c("try_job", "design_by_scheme", "SEAFAR_INIT")) %dopar% {
      try_job(job)
    }
  )

  stopCluster(cl)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  cat(sprintf("Parallel run finished in %s mins across %d worker(s).\n", elapsed, N_CORES))

  is_bad <- vapply(results_list, function(x) is.null(x) || inherits(x, "error"), logical(1))
  n_failed <- sum(is_bad)
  if (n_failed > 0) {
    cat(sprintf("%d / %d jobs failed and were excluded (see seafar_errors.log for the ones caught inside try_job()).\n",
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
  write.csv(results_df, "Seafar_results.csv", row.names = FALSE)

  design_group_vars <- c("Scheme", "Folder", "n_variables", "s_size", "p_sparse", "n_components", "VAFx")
  design_group_vars <- intersect(design_group_vars, names(results_df))

  final_summary <- results_df %>%
    group_by(across(all_of(design_group_vars))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  write.csv(final_summary, "Seafar_summary.csv", row.names = FALSE)

  # Stability against trueW; add "_vsP" metric names too if you also want
  # between-replication stability of the vs-trueP comparison.
  stability_summary <- summarise_stability(results_df, group_vars = design_group_vars)
  write.csv(stability_summary, "Seafar_stability.csv", row.names = FALSE)

  stability_summary_vsP <- summarise_stability(
    results_df, group_vars = design_group_vars,
    metric_vars = c("Recovery_Rate_vsP", "Precision_vsP", "Recall_vsP", "F1_vsP", "FPR_vsP", "FNR_vsP")
  )
  write.csv(stability_summary_vsP, "Seafar_stability_vsP.csv", row.names = FALSE)

  cat("Process finished successfully!\n")
}
