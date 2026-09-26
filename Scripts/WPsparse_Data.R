# P and W sparse using foreach
# Rosember Guerra
# 15-09-2019

# install.packages("doParallel")
# install.packages("mvtnorm")
# install.packages("MASS")

rm(list = ls(all.names = TRUE))

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

dir.create("DATA-R-WP-Sparse", showWarnings = FALSE) # Directory to save the data

# setting the number of cores 
library(doParallel)

no_cores <- detectCores() - 1
c1 <- makePSOCKcluster(no_cores)
registerDoParallel(c1)

# sparse PCA data simulation #
VAFx = c(.80, .95,1)           # Proportion of explained variance
p_sparse = c(.7,.8,.9)      # Proportion of sparsity
n_components = c(2,3)         # Number of components
s_size = c(100,500)           # Sample size
n_variables = c(10,100,1000)     # Number of variables
n_replications = c(100)        # Number of repetitions

design_matrix <- expand.grid( n_variables=n_variables,s_size=s_size,p_sparse=p_sparse,
                              n_components=n_components,VAFx =VAFx)

design_matrix_replication <- design_matrix[rep(1:nrow(design_matrix), times = n_replications), ]

Info_simulation = list(n_data_sets = nrow(design_matrix_replication), n_replications  =n_replications,
                        design_matrix_replication = design_matrix_replication)
save(Info_simulation, file = "DATA-R-WP-Sparse/Info_simulation.RData")

# start simulating the data 

# NOTE on the two branches saved per replication:
#  - X / Xtrue: loadings P are sparse (disjoint support across components,
#    by construction of indxnonzero below); the true weight matrix W is the
#    minimum-norm solution of Xtrue %*% W = tmat, generally dense. This is
#    the "P sparse, W free" condition.
#  - X2 / Xtrue2: generated under the *symmetric* model X = X W2 W2', i.e.
#    weights and loadings are literally the same sparse matrix W2. This is
#    the "W = P, both sparse" condition (the XWW' case discussed with
#    Katrijn -- see project literature note).
# Previously only W2 (for X2) was saved; W (for X) was not, which is why
# Demo/TSPCA_Sim.R needed manual "use out$W2 instead of out$W" swapping. Both are
# now saved, with consistent names, so downstream scripts can pick the right
# branch by name instead of hand-editing.

results_sim1_data1 <- foreach(i=1:nrow(design_matrix_replication),
                              .options.RNG = 2018,
                              .packages = c("MASS","mvtnorm"),
                              .combine=rbind,
                              .errorhandling = "pass")%dopar%{
                              tryCatch({

                                # Define the file path for this iteration
                                file_path <- paste0("DATA-R-WP-Sparse/WPsparse", i, ".RData")

                                # MODIFIED CHECK: Only skip if file exists AND is NOT 0 bytes
                                if (file.exists(file_path) && file.info(file_path)$size > 0) {
                                  return(NULL) # Move to the next iteration safely
                                }

                                # set the specific values of the parameters
                                R = design_matrix_replication$n_components[i]
                                I =  design_matrix_replication$s_size[i]
                                J = design_matrix_replication$n_variables[i]
                                vafx = design_matrix_replication$VAFx[i]
                                propsparse = design_matrix_replication$p_sparse[i]

                                n_zeros = floor(J*propsparse)
                                n_non_zero = J-n_zeros

                                # indxnonzero draws n_non_zero*R DISTINCT indices out of 1:J and
                                # reshapes them into R disjoint columns, so no variable ever loads on
                                # more than one component. That only works if there are enough
                                # variables to go around; fail with a clear message (caught and logged
                                # below) rather than the cryptic "cannot take a sample larger than the
                                # population" error base R would otherwise throw.
                                if (n_non_zero * R > J) {
                                  stop(sprintf(
                                    "n_non_zero (%d) * R (%d) = %d exceeds J (%d): cannot draw disjoint nonzero index sets for this design cell (p_sparse=%.2f). Lower R, raise p_sparse, or allow overlapping support across components.",
                                    n_non_zero, R, n_non_zero * R, J, propsparse))
                                }

                                # generating the data
                                X = mvrnorm(n =I,mu=rep(0,J),Sigma = diag(1,J))

                                # SVD on Xinit
                                svd1 <- svd(X)

                                # P matrix is the right singular vectors
                                P <- svd1$v[,1:R]
                                indxnonzero = matrix(sample(1:J, size = (n_non_zero*R)),n_non_zero,R )
                                for (z in 1:R) {
                                  P[-indxnonzero[,z],z] = 0
                                }

                                # Normalizing the columns of P
                                W2  <- P %*% diag(1/sqrt(diag(t(P) %*% P)))
                                P <- P%*%diag(svd1$d[1:R])
                                tmat <- svd1$u[,1:R]
                                tmat2 = X%*%W2
                                Xtrue <- tmat %*% t(P)
                                Xtrue2 = tmat2%*%t(W2)

                                # True weights for the X / Xtrue branch (see note above): the
                                # minimum-norm W solving Xtrue %*% W = tmat exactly, same approach as
                                # Psparse_Data.R. Xtrue has rank R < J so this system is underdetermined;
                                # ginv() returns the minimum Frobenius-norm solution.
                                W <- MASS::ginv(Xtrue) %*% tmat

                                # Adding noice 1 %
                                SSqXtrue =  sum(Xtrue^2)                        # sum squares of the data set
                                EX = mvrnorm(I,mu= rep(0,J),Sigma = diag(1,J))  # EX = Error of X
                                SSqEX = sum(EX^2)                               # Sum squares fo the EX
                                fx = sqrt(SSqXtrue*(1-vafx)/(vafx * SSqEX))
                                Xnew = Xtrue + fx*EX                            # Data with noise

                                # Adding noice 2 %
                                SSqXtrue =  sum(Xtrue2^2)                        # sum squares of the data set
                                EX = mvrnorm(I,mu= rep(0,J),Sigma = diag(1,J))  # EX = Error of X
                                SSqEX = sum(EX^2)                               # Sum squares fo the EX
                                fx = sqrt(SSqXtrue*(1-vafx)/(vafx * SSqEX))
                                Xnew2 = Xtrue2 + fx*EX                            # Data with noise

                                # Saving the data and structure
                                out = list(X = Xnew, W = W, X2=Xnew2, W2 = W2 ,P = P, Z = tmat,Z2=tmat2 ,k = R,
                                           Propsparse = propsparse, nzeros =  n_zeros, teller = i)
                                save(out, file=file_path)

                                return(NULL)

                              }, error = function(e) {
                                # Log and move on instead of aborting the whole parallel run (this is
                                # how the previously-missing WPsparse72.RData, etc. should surface next
                                # time instead of silently disappearing).
                                msg <- sprintf("[WPsparse %d] FAILED: %s", i, conditionMessage(e))
                                cat(msg, "\n", file = "DATA-R-WP-Sparse/generation_errors.log", append = TRUE)
                                return(msg)
                              })
                              }

# stop Cluster
stopCluster(c1)