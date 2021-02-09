# Plot SVD
eigplot <- function(A, dmax = 100) {
  require(rARPACK)
  eigres <- eigs(A = A, k = dmax)
  print(plot(sort(abs(eigres$values), decreasing = T), main = "Singular values of adjacency matrix", ylab = "Value"))
}

# Plot all SVDs
eigplots <- function(Alist, dmax = 100) {
  m <- length(Alist)
  cols <- floor(sqrt(m))
  rows <- ceiling(sqrt(m))
  par(mfrow = c(rows, cols))
  sapply(Alist, eigplot, dmax = dmax)
}


# Plot SVD
svdplot <- function(A, dmax = 100) {
  require(rARPACK)
  svdres <- svds(A = A, k = dmax)
  print(plot(svdres$d, main = "Singular values of adjacency matrix", ylab = "Value"))
}

# Plot all SVDs
svdplots <- function(Alist, dmax = 100) {
  m <- length(Alist)
  cols <- floor(sqrt(m))
  rows <- ceiling(sqrt(m))
  par(mfrow = c(rows, cols))
  sapply(Alist, eigplot, dmax = dmax)
}
