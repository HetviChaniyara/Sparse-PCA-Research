# Block Generalized Power Method for Sparse PCA (l0 penalty)
# Katrijn van Deun, September 2026
#
# Implements the block-l0 variant of gPower (Journee, Nesterov, Richtarik &
# Sepulchre, 2010, "Generalized Power Method for Sparse Principal Component
# Analysis", JMLR 11, 517-553, Section 3.2), plus a direct-cardinality
# variant (gpower_block_cardinality()) so it can be dropped into the same
# benchmark loops as elasticnet::spca(..., para = phi) and
# CEC_PLS_SEM(..., phi = phi) elsewhere in this project (see
# SPCA_Functions.R, Demo/Elastic_net.R, Demo/SPCA.R).
#
# Problem solved by gpower_block_l0() (A = data, I x J, optionally
# centered/scaled; U's columns orthonormal, I x R "generalized score
# directions"; Z's columns unit-norm sparse loadings, J x R):
#
#   max_{U'U = I_R, ||z_j||=1}  sum_j  mu_j^2 (u_j' A z_j)^2 - gamma_j ||z_j||_0
#
# For fixed U, the inner maximization over each z_j is closed-form hard
# thresholding of the correlations A'u_j. For fixed Z, U is updated via the
# polar decomposition (SVD) of A Z diag(mu). Alternating the two increases
# the objective monotonically (Journee et al. 2010, Prop. 8), which is used
# as a correctness check in test_gpower_block_l0.R.
#
# gpower_block_l0() returns $weights/$loadings/$scores in the same shape as
# CEC_PLS_SEM() and elasticnet::spca(), so it drops into the existing
# benchmark loops (align_components(), evaluate_variable_selection(),
# compute_vaf()) unchanged.

#' Per-component gamma reference values for gpower_block_l0()
#'
#' Returns the squared singular values of X (after the same centering/scaling
#' gpower_block_l0() would apply), one per component. Multiply by a relative
#' sparsity rho_j in (0, 1) to get an absolute gamma_j that accounts for
#' later components naturally having smaller achievable correlations than
#' the first -- see Details in \code{\link{gpower_block_l0}}.
#'
#' @inheritParams gpower_block_l0
#' @return numeric vector of length R.
#' @export
gpower_gamma_reference <- function(X, R, center = TRUE, scale = FALSE) {
  A <- scale(as.matrix(X), center = center, scale = scale)
  svd(A, nu = 0, nv = 0)$d[seq_len(R)]^2
}

#' Block Generalized Power Method for Sparse PCA (l0 penalty)
#'
#' @param X data matrix, I (observations) x J (variables).
#' @param R number of components.
#' @param gamma sparsity threshold(s), in squared-correlation units (see
#'   Details): a single value recycled to all R components, or a length-R
#'   vector. Larger gamma gives sparser loadings; gamma = 0 recovers
#'   ordinary PCA (up to the mu-weighted block ambiguity, see Details).
#' @param mu nonnegative per-component block weight(s), a single value
#'   (recycled) or length-R vector. Only the *relative* weighting across
#'   components matters. Default 1 for all components -- the paper requires
#'   mu_1 >= ... >= mu_R > 0 to fully break rotational ambiguity between
#'   equally-weighted components, but equal mu (the default) still converges
#'   to a valid stationary point, as gpowerr's own default confirms.
#' @param center,scale passed to base::scale() before fitting.
#' @param max_iter maximum number of power-method iterations.
#' @param tol relative-change convergence tolerance on the objective value.
#' @param verbose print the objective value at every iteration.
#'
#' @details
#' gamma is an absolute threshold: index i is kept in component j's support
#' iff (mu_j * (A'u_j)_i)^2 > gamma_j. Do NOT reference every component's
#' gamma to the same global scale (e.g. the largest singular value): later
#' components explain less variance and so have much smaller achievable
#' correlations, and a gamma calibrated to component 1 will zero out
#' component 3 long before component 1 is meaningfully sparse. Instead use
#' \code{\link{gpower_gamma_reference}} to get one reference value per
#' component (that component's own squared singular value from an initial
#' SVD) and scale each by its own relative sparsity rho_j in (0, 1):
#' \code{gamma <- rho * gpower_gamma_reference(X, R, center, scale)}.
#'
#' Known limitation: because the support of each z_j is a discontinuous
#' (hard-thresholded) function of U, a coefficient sitting very close to the
#' sqrt(gamma_j) boundary can flip in and out of the support across a few
#' iterations, producing tiny non-monotonic wobbles in the objective late in
#' convergence (checked down to tol = 1e-10 on simulated data: wobbles stay
#' on the order of 1e-4 relative to the objective and always settle). This
#' is an inherent property of l0 (as opposed to l1/soft-threshold) block
#' coordinate ascent, not a sign of divergence; the default tol = 1e-4 stops
#' comfortably before it becomes visible.
#'
#' @return list with:
#'   \item{weights}{J x R matrix of unit-norm sparse loadings (Z)}
#'   \item{loadings}{J x R orthonormal reconstruction loadings (P), obtained
#'     as in Elastic_net.R: P = polar(X'X %*% weights)}
#'   \item{scores}{I x R matrix, X %*% weights}
#'   \item{exp_var}{proportion of variance of X accounted for by weights/loadings}
#'   \item{sparsity}{proportion of zero entries, per component}
#'   \item{iterations}{number of power-method iterations used}
#'   \item{objective}{objective value at each iteration (for diagnosing convergence)}
#' @references Journee, M., Nesterov, Y., Richtarik, P. and Sepulchre, R.
#'   (2010). Generalized Power Method for Sparse Principal Component
#'   Analysis. Journal of Machine Learning Research, 11, 517-553.
#' @export
gpower_block_l0 <- function(X, R, gamma, mu = 1, center = TRUE, scale = FALSE,
                             max_iter = 1000, tol = 1e-4, verbose = FALSE) {

  A <- scale(as.matrix(X), center = center, scale = scale)
  I <- nrow(A)
  J <- ncol(A)

  if (R < 1 || R > min(I, J)) stop("R must be between 1 and min(nrow(X), ncol(X)).")
  if (length(gamma) == 1) gamma <- rep(gamma, R)
  if (length(gamma) != R) stop("gamma must have length 1 or R.")
  if (any(gamma < 0)) stop("gamma must be non-negative.")
  if (length(mu) == 1) mu <- rep(mu, R)
  if (length(mu) != R) stop("mu must have length 1 or R.")
  if (any(mu <= 0)) stop("mu must be strictly positive.")

  # Deterministic warm start: leading left singular vectors of A, i.e. the
  # unpenalized solution -- avoids the random-init reproducibility hazard of
  # the original MATLAB/gpowerr initializer.
  U <- svd(A, nu = R, nv = 0)$u

  gamma_mat <- matrix(gamma, nrow = J, ncol = R, byrow = TRUE)
  obj_trace <- numeric(max_iter)
  iter <- 0
  converged <- FALSE

  repeat {
    iter <- iter + 1

    Y <- t(A) %*% U                    # J x R, Y[, j] = A' u_j
    Ymu <- sweep(Y, 2, mu, `*`)        # scale column j by mu_j
    thresholded <- pmax(Ymu^2 - gamma_mat, 0)
    obj_trace[iter] <- sum(thresholded)

    if (obj_trace[iter] == 0) {
      stop("gamma is too large: every component's support is empty. ",
           "Reduce gamma (try values well below max(svd(X)$d)^2).")
    }

    # X-update is the polar decomposition of A %*% Z %*% diag(mu), where Z's
    # columns are the actual unit-norm thresholded loadings -- using the raw
    # (unnormalized) correlations here instead would make the per-column
    # scale of the gradient depend on each pattern's incidental norm as well
    # as mu_j, which for unequal mu can produce a slightly non-monotone step
    # right at convergence (a subtlety the mu=1 case happens to be immune
    # to, since all columns are then scaled identically either way).
    grad <- matrix(0, I, R)
    for (j in seq_len(R)) {
      pattern_j <- thresholded[, j] > 0
      if (!any(pattern_j)) {
        grad[, j] <- U[, j]  # this component's support is empty; leave it in place
      } else {
        z_j <- Y[pattern_j, j] / sqrt(sum(Y[pattern_j, j]^2))
        grad[, j] <- mu[j] * (A[, pattern_j, drop = FALSE] %*% z_j)
      }
    }
    svd_grad <- svd(grad)
    U <- svd_grad$u %*% t(svd_grad$v)

    if (verbose) message("iter ", iter, ": objective = ", obj_trace[iter])

    if (iter > 1) {
      rel_change <- abs(obj_trace[iter] - obj_trace[iter - 1]) / abs(obj_trace[iter - 1])
      if (rel_change < tol) { converged <- TRUE; break }
    }
    if (iter >= max_iter) break
  }
  if (!converged) warning("gpower_block_l0 reached max_iter (", max_iter, ") without converging.")

  # Final extraction: the support comes from the mu-scaled correlations, but
  # the loading VALUES are the raw (unscaled) correlations restricted to that
  # support -- mu only ever weights how much each component's objective term
  # counts towards the joint sum, it does not rescale the loading itself.
  Y <- t(A) %*% U
  Ymu <- sweep(Y, 2, mu, `*`)
  pattern <- Ymu^2 > gamma_mat
  Z <- Y
  Z[!pattern] <- 0

  zero_col <- colSums(Z != 0) == 0
  if (any(zero_col)) {
    stop("Component(s) ", paste(which(zero_col), collapse = ", "),
         " collapsed to an all-zero loading; reduce gamma for those components.")
  }
  Z <- sweep(Z, 2, sqrt(colSums(Z^2)), `/`)

  scores <- A %*% Z

  # Reconstruction loadings via polar decomposition, matching the convention
  # already used for the elastic-net benchmark in Elastic_net.R.
  cross <- t(A) %*% A %*% Z
  svd_cross <- svd(cross)
  P <- svd_cross$u %*% t(svd_cross$v)

  # Inlined rather than calling SPCA_Functions.R::compute_vaf(), so this file
  # has no source() ordering dependency; identical formula (1 - SSE/SST).
  X_hat <- A %*% Z %*% t(P)
  exp_var <- 1 - sum((A - X_hat)^2) / sum(A^2)

  list(
    weights = Z,
    loadings = P,
    scores = scores,
    exp_var = exp_var,
    sparsity = colMeans(Z == 0),
    iterations = iter,
    objective = obj_trace[seq_len(iter)]
  )
}

#' Block Generalized Power Method for Sparse PCA (direct cardinality control)
#'
#' Same power-iteration structure as \code{\link{gpower_block_l0}} (U is
#' updated via the polar decomposition of A %*% Z %*% diag(mu)), but Z's
#' support is chosen by DIRECTLY selecting the largest-magnitude mu-weighted
#' correlations up to an exact nonzero budget phi, instead of a continuous
#' gamma threshold -- mirroring \code{apply_cardinality()} in
#' SPCA_Functions.R, right down to the same vocabulary: "per_component"
#' ranks each column separately and keeps exactly phi_j per column;
#' "total" ranks all J*R entries of the weight matrix together and keeps
#' only the phi largest overall, letting components compete freely for the
#' shared budget instead of each getting a fixed share.
#'
#' Because top-k selection is a closed-form maximizer of the fixed-U
#' subproblem (see Details), this hits the requested cardinality EXACTLY
#' every time -- no gamma calibration, no bisection search needed. Prefer
#' this function whenever what you actually want is an exact nonzero count
#' (which is what every other benchmark in this project specifies); reach
#' for the gamma-based \code{\link{gpower_block_l0}} instead when you want
#' continuous regularization-path behavior (e.g. a smooth sequence of fits
#' as gamma varies) rather than a fixed target count.
#'
#' @param X,mu,center,scale,max_iter,tol,verbose see \code{\link{gpower_block_l0}}.
#' @param R number of components.
#' @param phi target cardinality. For \code{cardinality_type = "per_component"}:
#'   a single value (recycled to every component) or a length-R vector, each
#'   in [1, ncol(X)] -- component j always ends up with exactly phi_j
#'   nonzero loadings. For \code{cardinality_type = "total"}: a single
#'   scalar in [R, R * ncol(X)] -- exactly phi nonzero loadings across the
#'   WHOLE J x R weight matrix in total, allocated freely across components
#'   by whichever correlations are largest.
#' @param cardinality_type "per_component" (default) or "total". As already
#'   noted for \code{apply_cardinality()}'s "total" mode in
#'   SPCA_Functions.R, ranking the whole matrix together CAN starve a
#'   component down to zero nonzero entries if its correlations are
#'   uniformly weaker than the others' -- this function stops with an
#'   informative error in that case rather than silently returning a
#'   degenerate component; try "per_component", a larger phi, or mu weights
#'   that favor the starved component.
#'
#' @details
#' For fixed U, maximizing sum_j mu_j^2 (u_j'Az_j)^2 over unit-norm z_j
#' subject to a cardinality budget has a closed-form solution: for any
#' fixed support S, the optimal z_j (restricted to S) is proportional to
#' the correlations Y_ij = (A'u_j)_i on S, so keeping entry (i,j) is worth
#' exactly (mu_j Y_ij)^2 to the objective, independent of which other
#' entries are kept. Under a per-column budget phi_j, that means: keep
#' column j's phi_j largest |mu_j Y_ij| values. Under a single TOTAL budget
#' phi shared across columns, entries are still mutually independent in
#' value, so the globally optimal allocation is exactly the phi largest
#' |mu_j Y_ij| values anywhere in the J x R matrix -- a separable-knapsack
#' argument (item (i,j)'s value never depends on what else is chosen, so
#' greedy-by-value is optimal), the same rule
#' \code{apply_cardinality(..., "total")} already applies to a weight
#' update in SPCA_Functions.R, here applied to the correlations instead.
#'
#' @return same shape as \code{\link{gpower_block_l0}}, plus \code{$phi} and
#'   \code{$cardinality_type} echoing what was requested.
#' @export
gpower_block_cardinality <- function(X, R, phi, mu = 1,
                                      cardinality_type = c("per_component", "total"),
                                      center = TRUE, scale = FALSE,
                                      max_iter = 1000, tol = 1e-4, verbose = FALSE) {

  cardinality_type <- match.arg(cardinality_type)
  A <- scale(as.matrix(X), center = center, scale = scale)
  I <- nrow(A)
  J <- ncol(A)

  if (R < 1 || R > min(I, J)) stop("R must be between 1 and min(nrow(X), ncol(X)).")
  if (length(mu) == 1) mu <- rep(mu, R)
  if (length(mu) != R) stop("mu must have length 1 or R.")
  if (any(mu <= 0)) stop("mu must be strictly positive.")

  if (cardinality_type == "per_component") {
    if (length(phi) == 1) phi <- rep(phi, R)
    if (length(phi) != R) stop("phi must have length 1 or R when cardinality_type = 'per_component'.")
    if (any(phi < 1 | phi > J)) stop("phi must be between 1 and ncol(X).")
  } else {
    if (length(phi) != 1) stop("phi must be a single number when cardinality_type = 'total'.")
    if (phi < R || phi > J * R) stop("phi must be between R and R * ncol(X) when cardinality_type = 'total'.")
  }

  # Deterministic warm start, as in gpower_block_l0() -- see the comment
  # there for why this matters beyond mere reproducibility.
  U <- svd(A, nu = R, nv = 0)$u

  select_pattern <- function(Ymu2) {
    pattern <- matrix(FALSE, J, R)
    if (cardinality_type == "per_component") {
      for (j in seq_len(R)) {
        keep <- order(Ymu2[, j], decreasing = TRUE)[seq_len(phi[j])]
        pattern[keep, j] <- TRUE
      }
    } else {
      keep <- order(Ymu2, decreasing = TRUE)[seq_len(phi)]
      pattern[keep] <- TRUE
    }
    pattern
  }

  check_no_empty_component <- function(pattern) {
    empty <- colSums(pattern) == 0
    if (any(empty)) {
      stop("Component(s) ", paste(which(empty), collapse = ", "),
           " received zero nonzero entries under this ",
           if (cardinality_type == "total") "total budget (phi = " else "phi (",
           paste(phi, collapse = ", "), "). ",
           if (cardinality_type == "total")
             "Try cardinality_type = 'per_component', a larger phi, or mu weights that favor the starved component(s)."
           else "Increase phi for the affected component(s).")
    }
  }

  obj_trace <- numeric(max_iter)
  iter <- 0
  converged <- FALSE

  repeat {
    iter <- iter + 1

    Y <- t(A) %*% U
    Ymu <- sweep(Y, 2, mu, `*`)
    Ymu2 <- Ymu^2
    pattern <- select_pattern(Ymu2)
    check_no_empty_component(pattern)
    obj_trace[iter] <- sum(Ymu2[pattern])

    grad <- matrix(0, I, R)
    for (j in seq_len(R)) {
      pattern_j <- pattern[, j]
      z_j <- Y[pattern_j, j] / sqrt(sum(Y[pattern_j, j]^2))
      grad[, j] <- mu[j] * (A[, pattern_j, drop = FALSE] %*% z_j)
    }
    svd_grad <- svd(grad)
    U <- svd_grad$u %*% t(svd_grad$v)

    if (verbose) message("iter ", iter, ": objective = ", obj_trace[iter])

    if (iter > 1) {
      rel_change <- abs(obj_trace[iter] - obj_trace[iter - 1]) / abs(obj_trace[iter - 1])
      if (rel_change < tol) { converged <- TRUE; break }
    }
    if (iter >= max_iter) break
  }
  if (!converged) warning("gpower_block_cardinality reached max_iter (", max_iter, ") without converging.")

  # Final extraction using the same selection rule as every iteration.
  Y <- t(A) %*% U
  Ymu2 <- sweep(Y, 2, mu, `*`)^2
  pattern <- select_pattern(Ymu2)
  check_no_empty_component(pattern)

  Z <- Y
  Z[!pattern] <- 0
  Z <- sweep(Z, 2, sqrt(colSums(Z^2)), `/`)

  scores <- A %*% Z

  cross <- t(A) %*% A %*% Z
  svd_cross <- svd(cross)
  P <- svd_cross$u %*% t(svd_cross$v)

  X_hat <- A %*% Z %*% t(P)
  exp_var <- 1 - sum((A - X_hat)^2) / sum(A^2)

  list(
    weights = Z,
    loadings = P,
    scores = scores,
    exp_var = exp_var,
    sparsity = colMeans(Z == 0),
    iterations = iter,
    objective = obj_trace[seq_len(iter)],
    phi = phi,
    cardinality_type = cardinality_type
  )
}
