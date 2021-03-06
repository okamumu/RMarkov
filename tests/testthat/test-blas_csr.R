library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for blas of CSR")

# this test requires src/test_blascsr.cpp

create.sparse <- function(m, n, p, min=-10, max=10) {
  np <- ceiling(m * n * p)
  A <- matrix(runif(m*n, min=min, max=max), m, n)
  i <- sample(x=1:m, size=np, replace=TRUE)
  j <- sample(x=1:n, size=np, replace=TRUE)
  for (u in 1:np) {
    A[i[u],j[u]] <- 0
  }
  A
}

dense_to_sparse <- function(x) {
  as(Matrix(x), "RsparseMatrix")
}

test_that("blas csr dgemvN", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(20, min = -10, max = 10)
  y <- runif(10, min = -10, max = 10)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  spA <- dense_to_sparse(A)
  z <- alpha * as.vector(A %*% x) + beta * y
  C_dgemvN_csr(alpha, spA, x, beta, y)
  expect_equal(y, z)
})

test_that("blas csr dgemvT", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(10, min = -10, max = 10)
  y <- runif(20, min = -10, max = 10)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  spA <- dense_to_sparse(A)
  z <- alpha * as.vector(t(A) %*% x) + beta * y
  C_dgemvT_csr(alpha, spA, x, beta, y)
  expect_equal(y, z)
})

test_that("blas csr dgerN", {
  alpha <- runif(1)
  m <- 20
  n <- 10
  x <- runif(m, min = -10, max = 10)
  y <- runif(n, min = -10, max = 10)
  X <- matrix(x, ncol=1) %*% matrix(y, nrow=1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(m, n, p)
  X[A == 0] <- 0
  spA <- dense_to_sparse(A)
  z <- alpha * X + A
  spZ <- dense_to_sparse(z)
  C_dgerN_csr(alpha, x, y, spA)
  expect_equal(spA, spZ)
})

test_that("blas csr dgerT", {
  alpha <- runif(1)
  m <- 20
  n <- 10
  x <- runif(m, min = -10, max = 10)
  y <- runif(n, min = -10, max = 10)
  X <- matrix(x, ncol=1) %*% matrix(y, nrow=1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(m, n, p)
  X[A == 0] <- 0
  spA <- dense_to_sparse(A)
  z <- alpha * X + A
  spZ <- dense_to_sparse(z)
  C_dgerT_csr(alpha, y, x, spA)
  expect_equal(spA, spZ)
})

test_that("blas csr dgemmNN", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(m, k, p)
  spA <- dense_to_sparse(A)
  B <- matrix(runif(k*n, min = -10, max = 10), k, n)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (A %*% B) + beta * C
  C_dgemmNN_csr(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmTN", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(k, m, p)
  spA <- dense_to_sparse(A)
  B <- matrix(runif(k*n, min = -10, max = 10), k, n)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (t(A) %*% B) + beta * C
  C_dgemmTN_csr(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmNT", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(m, k, p)
  spA <- dense_to_sparse(A)
  B <- matrix(runif(k*n, min = -10, max = 10), n, k)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (A %*% t(B)) + beta * C
  C_dgemmNT_csr(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmTT", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(k, m, p)
  spA <- dense_to_sparse(A)
  B <- matrix(runif(k*n, min = -10, max = 10), n, k)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (t(A) %*% t(B)) + beta * C
  C_dgemmTT_csr(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmNN2", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(m*k, min = -10, max = 10), m, k)
  p <- runif(1, min=0.1, max=0.9)
  B <- create.sparse(k, n, p)
  spB <- dense_to_sparse(B)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (A %*% B) + beta * C
  C_dgemmNN2_csr(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmTN2", {
  m <- 20
  n <- 10
  k <- 5
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(m*k, min = -10, max = 10), k, m)
  p <- runif(1, min=0.1, max=0.9)
  B <- create.sparse(k, n, p)
  spB <- dense_to_sparse(B)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (t(A) %*% B) + beta * C
  C_dgemmTN2_csr(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmNT2", {
  m <- 2
  n <- 1
  k <- 5
  alpha <- 1 #runif(1)
  beta <- 0 #runif(1)
  A <- matrix(runif(m*k, min = -10, max = 10), m, k)
  p <- runif(1, min=0.1, max=0.9)
  B <- create.sparse(n, k, p)
  spB <- dense_to_sparse(B)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (A %*% t(B)) + beta * C
  C_dgemmNT2_csr(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csr dgemmTT2", {
  m <- 2
  n <- 1
  k <- 5
  alpha <- 1 #runif(1)
  beta <- 0 #runif(1)
  A <- matrix(runif(m*k, min = -10, max = 10), k, m)
  p <- runif(1, min=0.1, max=0.9)
  B <- create.sparse(n, k, p)
  spB <- dense_to_sparse(B)
  C <- matrix(runif(m*n, min = -10, max = 10), m, n)
  z <- alpha * (t(A) %*% t(B)) + beta * C
  C_dgemmTT2_csr(alpha, A, spB, beta, C)
  expect_equal(C, z)
})
