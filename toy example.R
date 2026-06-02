library(elasticnet)

X1 <- matrix(c(
  1.0,  1.0,  0.9,  0.0,  0.0,  0.0,
  -1.0, -1.0, -0.9,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,
  0.0,  0.0,  0.0, -1.0, -1.0, -1.0,
  0.5,  0.5,  0.5,  0.5,  0.5,  0.5
), nrow = 5, byrow = TRUE)

colnames(X1) <- paste0("V", 1:6)
rownames(X1) <- paste0("Obs", 1:5)

X2 <- matrix(c(
  1.0,  0.7,  1.2,  0.0,  0.0,  0.0,  
  -1.0, -0.7, -1.2,  0.0,  0.0,  0.0,
  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,
  0.0,  0.0,  0.0, -1.0, -1.0, -1.0,
  0.5,  0.5,  0.5,  0.5,  0.5,  0.5
), nrow = 5, byrow = TRUE)

colnames(X2) <- paste0("V", 1:6)
rownames(X2) <- paste0("Obs", 1:5)

# Classical PCA
# col's were flipped to match vars 1-3 to component 1
pca_model1 <- prcomp(X1, center = FALSE, scale. = FALSE) 
W_PCA <- pca_model1$rotation[, c(2, 1)]
colnames(W_PCA) <- c("Comp1", "Comp2")

# Zou's Sparse PCA 
spca_model1 <- arrayspc(X1, K = 2, para = c(2, 2), trace = FALSE)
spca_model2 <- arrayspc(X2, K = 2, para = c(2, 2), trace = FALSE)

W_sPCA_X1 <- spca_model1$loadings[, c(2, 1)]
W_sPCA_X2 <- spca_model2$loadings[, c(2, 1)]
colnames(W_sPCA_X1) <- colnames(W_sPCA_X2) <- c("Comp1", "Comp2")
rownames(W_sPCA_X1) <- rownames(W_sPCA_X2) <- colnames(X1)


run_uslpca <- function(X, K = 2, total_cardinality = 4, max_iter = 100, tol = 1e-7) {
  n <- nrow(X)
  p <- ncol(X)

  S <- (t(X) %*% X) / n 

  set.seed(123)
  A <- matrix(rnorm(p * K), p, K)
  A <- apply(A, 2, function(x) x / sqrt(sum(x^2)))
  
  ls_n_old <- 1 - (sum(A^2) / sum(diag(S)))
  
  for(iter in 1:max_iter) {
    ASA <- t(A) %*% S %*% A
    evd_out <- eigen(ASA, symmetric = TRUE)
    
    L <- evd_out$vectors

    Lambda_sq <- pmax(evd_out$values, 1e-10) 
    Lambda_inv <- 1 / sqrt(Lambda_sq)

    B <- S %*% A %*% L %*% (Lambda_inv * t(L))
    
    b_squared <- B^2
    ranks <- rank(-b_squared, ties.method = "first") 
    
    A_new <- B
    A_new[ranks > total_cardinality] <- 0
    
    ls_n_new <- 1 - (sum(A_new^2) / sum(diag(S)))
    
    if (abs(ls_n_old - ls_n_new) <= tol) {
      A <- A_new
      break
    }
    
    A <- A_new
    ls_n_old <- ls_n_new
  }
  
  return(A)
}

# Run USLPCA
P_USLPCA_X1 <- run_uslpca(X1, K = 2, total_cardinality = 6)[, c(2, 1)]
P_USLPCA_X2 <- run_uslpca(X2, K = 2, total_cardinality = 6)[, c(2, 1)]
colnames(P_USLPCA_X1) <- colnames(P_USLPCA_X2) <- c("Comp1", "Comp2")
rownames(P_USLPCA_X1) <- rownames(P_USLPCA_X2) <- colnames(X1)

# Outputs
cat("\n[1] Classical PCA:\n")
print(round(W_PCA, 2))

cat("\n[2] Zou'S sPCA X1:\n")
print(round(W_sPCA_X1, 2))

cat("\n[3] Zou'S sPCA ON X2:\n")
print(round(W_sPCA_X2, 2))

cat("\n[5] USLPCA X1:\n")
print(round(P_USLPCA_X1, 2))

cat("\n[6] USLPCA X2:\n")
print(round(P_USLPCA_X2, 2))