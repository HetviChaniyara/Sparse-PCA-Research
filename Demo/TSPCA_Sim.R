# TSPCA Research
# Updated Version March 2026 Hetvi Chaniyara
# Runs TSPCA for various data folders
# Incorporates the changes of Katrijn's December 2025 version
#
# Revised September 2026 (first pass): runs automatically over all
# data-generation schemes (W-sparse, P-sparse, WP-sparse) without
# hand-editing which `out$` field to read -- see the `conditions` table
# below. Also fixed: write.csv() previously received the raw list of
# per-dataset data.frames instead of the row-bound table; `if (prefix==1)`
# was always FALSE (prefix is a string) so recovery was evaluated against
# the wrong ground truth; folders/files that don't exist yet are now
# skipped with a clear message instead of crashing the whole run.
#
# Revised September 2026 (second pass): two further additions.
#  1. Stability of zero/nonzero status. Each replication in a design cell
#     draws its OWN random true support (see e.g. WPsparse_Data.R), so
#     "does item j stay zero/nonzero across replications" -- the approach
#     used for the fixed real Big5 data in Illustration_Big5.qmd, where the
#     same 240 items are fixed across resamples -- isn't a meaningful
#     question here: item j's true status is itself randomly redrawn every
#     replication. What IS meaningful, and computed here, is the
#     BETWEEN-REPLICATION VARIABILITY of the selection-quality metrics
#     within a design cell: if TSPCA reliably achieves similar precision/
#     recall/FPR/FNR every time it is run on a fresh dataset from that same
#     design, the between-replication SD/IQR of those metrics is small
#     (stable selection behaviour); if performance swings widely from
#     dataset to dataset despite an identical design, the SD/IQR is large
#     (unstable selection behaviour). See TSPCA_Functions.R::
#     evaluate_variable_selection() (now also returns FP/FN rates) and
#     ::summarise_stability() (does the grouped mean/SD/IQR aggregation).
#     Written to TSPCA_stability.csv, one row per design cell.
#  2. Parallel execution. The main loop now runs under doParallel/foreach,
#     matching the convention already used in Wsparse_Data.R /
#     WPsparse_Data.R, instead of looping sequentially over every dataset in
#     every scheme. See the "parallel execution" section below and the
#     accompanying note to Katrijn on tuning N_CORES / chunking for a
#     10,800-dataset-per-scheme grid.

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

source("../Scripts/TSPCA_Functions.R")
library(dplyr)
library(doParallel)
library(foreach)

# ---- configuration --------------------------------------------------------
# "total": cardinality budget is shared across the whole JxR weight matrix
#   (the original behaviour -- can allocate unevenly across components).
# "per_component": each component gets its own budget, derived the same way
#   SPCA_Sim.R already derives its per-column oracle target. Set this to
#   "per_component" for a fairer, crash-safer comparison; see
#   TSPCA_Functions.R::apply_cardinality for the mechanics.
CARDINALITY_TYPE <- "total"
CONSTRAINED      <- 0     # 0 = penalized (W minus P), 1 = equality-constrained (W = P)
N_MULTISTART     <- 10
MAX_ITER         <- 100

# Number of parallel worker processes. detectCores() - 1 leaves one core
# free for the OS/RStudio, matching the convention already used in
# Wsparse_Data.R / WPsparse_Data.R. Override with a smaller number if you
# want to keep using the machine for other work while this runs, or a
# larger one (up to detectCores()) on a dedicated run.
N_CORES <- max(1, parallel::detectCores() - 1)

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

# Reads the design matrix that the data-generation scripts save alongside the
# datasets in Info_simulation.RData.
load_design <- function(folder) {
  env <- new.env()
  load(file.path(folder, "Info_simulation.RData"), envir = env)
  env$Info_simulation$design_matrix_replication
}

# ---- build a flat job list across ALL schemes ------------------------------
# Flattening every (scheme, dataset-index) pair into one job list, rather
# than parallelising each scheme's loop separately one after another, keeps
# every core busy across the whole run: schemes with fewer datasets, or that
# happen to finish early, don't leave cores idle while another scheme is
# still going. Existence checks (folder, then each file) happen here,
# sequentially, before any worker is spun up -- cheap relative to a single
# TSPCA fit, and lets us report missing-file counts up front instead
# of discovering them mid-run inside a worker.
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
# Each job is one dataset's full N_MULTISTART-restart TSPCA fit --
# a natural unit of work (not too fine-grained, so per-task overhead stays
# negligible relative to the ~10 fits it does; not so coarse that load
# balancing suffers). Do NOT also parallelise the N_MULTISTART loop inside a
# job: PSOCK workers spawning their own sub-clusters is fragile and
# unnecessary here, since with thousands of jobs >> N_CORES, per-job
# parallelism alone already keeps every core saturated.
#
# clusterEvalQ() sources TSPCA_Functions.R and loads dplyr ONCE per worker at
# startup, rather than re-sourcing it on every one of the (potentially many
# thousands of) jobs -- PSOCK workers start with a fresh, empty environment
# and don't inherit anything from this script, so this step is required (a
# job that calls TSPCA() without it will fail with "could not find
# function").
#
# .errorhandling = "pass" (matching the convention already used in
# Wsparse_Data.R / WPsparse_Data.R) means one failing job returns its error
# instead of aborting the whole run; failures are also caught explicitly
# inside try_job() and logged to simulation_errors.log so they're easy to
# find afterwards, and filtered out of the final results below.
if (length(jobs) > 0) {

  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    source("../Scripts/TSPCA_Functions.R")
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
        round((1 - design$p_sparse[job$i]) * J) # recycled to every column by TSPCA
      }
      rho <- sum(X^2) / R

      # run with multistart
      best_res <- NULL; best_loss <- Inf
      for (m in 1:N_MULTISTART) {
        set.seed(100 + m)
        res <- TSPCA(X, INIT = NULL, R, 1e-8, phi, rho, constrained = CONSTRAINED,
                            MaxIter = MAX_ITER, cardinality_type = CARDINALITY_TYPE)
        if (res$Residual < best_loss) { best_loss <- res$Residual; best_res <- res }
      }

      # metrics calculation -- always against this condition's own ground
      # truth, no manual "Change to W2 if WPsparse" swapping needed
      W_aligned <- align_components(best_res$weights, trueW)
      P_aligned <- align_components(best_res$loadings, trueP)
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
        Loss = best_res$Residual,
        FEV = compute_vaf(X, best_res$weights, best_res$loadings),
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
        Iterations = best_res$n_iterations,
        MSE_W_P = bvm_W_P$mse
      )
    }, error = function(e) {
      msg <- sprintf("[%s / %s%d] FAILED: %s", job$scheme, job$prefix, job$i, conditionMessage(e))
      cat(msg, "\n", file = "simulation_errors.log", append = TRUE)
      NULL
    })
  }

  # foreach's automatic export turns out to be transitive: it doesn't just
  # scan the free variables of the expression literally inside the %dopar%
  # block (`try_job`, and the iteration variable `job`), it also walks into
  # try_job()'s own body (via codetools::findGlobals) to pick up the globals
  # IT reads -- design_by_scheme, N_MULTISTART, CARDINALITY_TYPE,
  # CONSTRAINED, MAX_ITER -- and ships exactly those to each worker (tested
  # directly: an unrelated global left out of that chain is NOT visible on
  # the worker, so this isn't "export all of .GlobalEnv" either). That
  # matters because try_job() is defined at this script's top level
  # (environment = .GlobalEnv on the MASTER process), and R's serialization
  # treats .GlobalEnv as a well-known special environment that gets remapped
  # to each WORKER's own (separate, near-empty) .GlobalEnv on arrival --
  # so anything NOT auto-exported would fail on the workers with an "object
  # not found" error despite running fine sequentially.
  #
  # The explicit .export list below is redundant with that auto-detection
  # (confirmed empirically) and is kept only as a belt-and-suspenders
  # safeguard in case a future refactor breaks the transitive chain (e.g.
  # try_job stops being a plain top-level function). Because both mechanisms
  # agree, foreach emits a harmless "already exporting variable(s)" warning
  # -- suppressed here since it doesn't indicate a problem.
  results_list <- suppressWarnings(
    foreach(job = jobs, .errorhandling = "pass",
            .export = c("try_job", "design_by_scheme", "N_MULTISTART",
                        "CARDINALITY_TYPE", "CONSTRAINED", "MAX_ITER")) %dopar% {
      try_job(job)
    }
  )

  stopCluster(cl)

  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  cat(sprintf("Parallel run finished in %s mins across %d worker(s).\n", elapsed, N_CORES))

  # Drop NULLs (jobs caught by try_job's own tryCatch) and any raw error
  # objects (jobs that failed outside try_job's tryCatch, e.g. a worker
  # crash -- .errorhandling = "pass" returns the condition object itself in
  # that case) before binding into one table.
  is_bad <- vapply(results_list, function(x) is.null(x) || inherits(x, "error"), logical(1))
  n_failed <- sum(is_bad)
  if (n_failed > 0) {
    cat(sprintf("%d / %d jobs failed and were excluded (see simulation_errors.log for the ones caught inside try_job()).\n",
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
  write.csv(results_df, "TSPCA_results.csv", row.names = FALSE)

  design_group_vars <- c("Scheme", "Folder", "n_variables", "s_size", "p_sparse", "n_components", "VAFx")
  design_group_vars <- intersect(design_group_vars, names(results_df))

  final_summary <- results_df %>%
    group_by(across(all_of(design_group_vars))) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  write.csv(final_summary, "TSPCA_summary.csv", row.names = FALSE)

  # Stability of zero/nonzero selection status across replications within
  # each design cell -- see the file header and
  # TSPCA_Functions.R::summarise_stability() for what this does and does not
  # capture in a simulation-study (as opposed to fixed-real-data
  # resampling) setting.
  stability_summary <- summarise_stability(results_df, group_vars = design_group_vars)
  write.csv(stability_summary, "TSPCA_stability.csv", row.names = FALSE)

  cat("Process finished successfully!\n")
}
