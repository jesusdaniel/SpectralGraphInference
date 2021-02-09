#'
#' Function to sample MASE statistic from the COSIE model. A sample of 2B graphs
#' with E[A] = VRV^T is obtained, and for every pair of them, the distance between their
#' score matrices obtained by MASE is calculated.
#'
#' @param V common invariant subspace of COSIE
#' @param R score matrix of COSIE
#' @param P_0 optional average matrix in the model. Default is 0.
#' @param B number of test statistics samples
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NULL, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augment logial. When true, the diagonal augmentation method for ASE is performed.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#'
#' @return A vector containing a sample of test statistics
#'
#' @references
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
mase.test.statistic.distr <- function(V, R, P0 = NULL, B = 1000,
                                      d = NA, d_vec = NA, scaled.ASE = FALSE, diag.augment = TRUE) {
  n = nrow(V)
  d = ncol(V)

  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1 <- sample_from_P(P1)
    Ab2 <- sample_from_P(P1)
    mase.embedding <- mase(list(Ab1-P0, Ab2-P0), d = d, d_vec = d_vec,
                           scaled.ASE = scaled.ASE, diag.augment = diag.augment)
    return(norm(mase.embedding$R[[1]]- mase.embedding$R[[2]], "F"))
  })
  return(bootstrap)
}

#'
#' Function to sample MASE vertex statistics from the COSIE model. A sample of 2B graphs
#' with E[A] = VRV^T is obtained, and for every pair of them, the distance between their
#' score matrices obtained by MASE is calculated.
#'
#' @param V common invariant subspace of COSIE
#' @param R score matrix of COSIE
#' @param P_0 optional average matrix in the model. Default is 0.
#' @param B number of test statistics samples
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NULL, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augment logial. When true, the diagonal augmentation method for ASE is performed.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#'
#' @return A vector containing a sample of test statistics
#'
#' @references
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
mase.test.statistic.vertex.distr <- function(V, R, P0 = NULL, B = 1000,
                                      d = NA, d_vec = NA, scaled.ASE = FALSE, diag.augment = TRUE) {
  n = nrow(V)
  d = ncol(V)

  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- lapply(1:B, function(i) {
    Ab1 <- sample_from_P(P1)
    Ab2 <- sample_from_P(P1)
    mase.embedding <- mase(list(Ab1-P0, Ab2-P0), d = d, d_vec = d_vec,
                           scaled.ASE = scaled.ASE, diag.augment = diag.augment)
    Xhats <- latent_positions(mase.embedding)
    Xdiff <- Xhats[[1]] - Xhats[[2]]
    norms <- apply(Xdiff, 1, function(x) sqrt(sum(x^2)))
    return(norms)
  })
  vertex.boostrap <- lapply(1:n, function(i) sapply(bootstrap, function(x) x[i]))
  return(vertex.boostrap)
}

latent_positions <- function(mase.object) {
  K <- dim(mase.object$V)[2]
  Xhats <- list()
  i <- 1
  for(R in mase.object$R) {
    eig = eigen(R)
    Rsqrt <- eig$vectors %*% diag(sign(eig$values) * sqrt(abs(eig$values)), K) %*%
      t(eig$vectors)
    Xhats[[i]] <- mase.object$V %*% Rsqrt
    i = i+1
  }
  return(Xhats)
}
#'
#' Given a set of test statistics and a sample null distribution, calculate p-values
p.values <- function(test.statistics, null.distribution) {
  return(sapply(test.statistics, function(test)
    (sum(test <= null.distribution) + 0.5)/ length(null.distribution)))
}


#'
#' Function to perform a parametric bootstrap to test the hypothesis
#' that two score matrices in the COSIE model are not different.
#'
#' @param V common invariant subspace of COSIE
#' @param R1 score matrix of COSIE for the first graph
#' @param R2 score matrix of COSIE for the first graph
#' @param B number of test statistics samples
#' @param P_0 optional average matrix in the model. Default is 0.
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NULL, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augment logical. When true, the diagonal augmentation method for ASE is performed.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#'
#' @return A list containing a matrix V of size n x d, with the
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#'
#' @references
#' Tang, Minh, et al. "A semiparametric two-sample hypothesis testing problem for random graphs."
#' Journal of Computational and Graphical Statistics 26.2 (2017): 344-354.
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
test_equality_COSIE <- function(V, R1, R2, B = 1000, P0 = NULL,
                                d = NA, d_vec = NA, scaled.ASE = FALSE, diag.augment = TRUE) {
  if(is.null(P0)) P0 = Matrix(0, n, n)
  null_1 <- mase.test.statistic.distr(V, R1, P0 = P0, B = B,
                                      d = d, d_vec = d_vec, scaled.ASE = scaled.ASE, diag.augment = diag.augment)
  null_2 <- mase.test.statistic.distr(V, R2, P0 = P0, B = B,
                                      d = d, d_vec = d_vec, scaled.ASE = scaled.ASE, diag.augment = diag.augment)
  test <- norm(R1 - R2, type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2,
              test.statistic = test))
}



