library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for blas of CSC")

# this test requires src/test_blascsc.cpp

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
  as(Matrix(x), "CsparseMatrix")
}

test_that("blas csc dgemvN", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(20, min = -10, max = 10)
  y <- runif(10, min = -10, max = 10)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  spA <- dense_to_sparse(A)
  z <- alpha * as.vector(A %*% x) + beta * y
  C_dgemvN_csc(alpha, spA, x, beta, y)
  expect_equal(y, z)
})

test_that("blas csc dgemvT", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(10, min = -10, max = 10)
  y <- runif(20, min = -10, max = 10)
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  spA <- dense_to_sparse(A)
  z <- alpha * as.vector(t(A) %*% x) + beta * y
  C_dgemvT_csc(alpha, spA, x, beta, y)
  expect_equal(y, z)
})

test_that("blas csc dgerN", {
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
  C_dgerN_csc(alpha, x, y, spA)
  expect_equal(spA, spZ)
})

test_that("blas csc dgerT", {
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
  C_dgerT_csc(alpha, y, x, spA)
  expect_equal(spA, spZ)
})

test_that("blas csc dgemmNN", {
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
  C_dgemmNN_csc(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmTN", {
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
  C_dgemmTN_csc(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmNT", {
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
  C_dgemmNT_csc(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmTT", {
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
  C_dgemmTT_csc(alpha, spA, B, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmNN2", {
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
  C_dgemmNN2_csc(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmTN2", {
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
  C_dgemmTN2_csc(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmNT2", {
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
  C_dgemmNT2_csc(alpha, A, spB, beta, C)
  expect_equal(C, z)
})

test_that("blas csc dgemmTT2", {
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
  C_dgemmTT2_csc(alpha, A, spB, beta, C)
  expect_equal(C, z)
})
