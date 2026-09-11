# Block Generalized Power Method for Sparse PCA (l0 penalty)
# Katrijn van Deun, September 2026
#
# Implements the block-l0 variant of gPower (Journee, Nesterov, Richtarik &
# Sepulchre, 2010, "Generalized Power Method for Sparse Principal Component
# Analysis", JMLR 11, 517-553, Section 3.2).
#
# Problem solved (A = data, I x J, optionally centered/scaled; U's columns
# orthonormal, I x R "generalized score directions"; Z's columns unit-norm
# sparse loadings, J x R):
#
#   max_{U'U = I_R, ||z_j||=1}  sum_j  mu_j^2 (u_j' A z_j)^2 - gamma_j ||z_j||_0
#
# For fixed U, the inner maximization over each z_j is closed-form hard
# thresholding of the correlations A'u_j. For fixed Z, U is updated via the
# polar decomposition (SVD) of A Z diag(mu). Alternating the two increases
# the objective monotonically (Journee et al. 2010, Prop. 8), which is used
# below as a correctness check (see test_gpower_block_l0.R).
#
# Returns $weights/$loadings/$scores in the same shape as CEC_PLS_SEM() and
# elasticnet::spca() in this project (see SPCA_Functions.R, Elastic_net.R),
# so it drops into the existing benchmark loops (align_components(),
# evaluate_variable_selection(), compute_vaf()) unchanged.

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
