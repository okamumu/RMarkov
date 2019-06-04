library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for blas2 (dense matrix)")

# this test requires src/test_blas2.cpp

test_that("blas2 dgemvN", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(2, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  z <- alpha * as.vector(A %*% x) + beta * y
  C_dgemvN(alpha, A, x, beta, y)
  expect_equal(y, z)
})

test_that("blas2 dgemvT", {
  alpha <- runif(1)
  beta <- runif(1)
  x <- runif(3, min = -10, max = 10)
  y <- runif(2, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  z <- alpha * as.vector(t(A) %*% x) + beta * y
  C_dgemvT(alpha, A, x, beta, y)
  expect_equal(y, z)
})

test_that("blas2 dgerN", {
  alpha <- runif(1)
  x <- runif(2, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 2, 3)
  z <- alpha * (matrix(x, ncol=1) %*% matrix(y, nrow=1)) + A
  C_dgerN(alpha, x, y, A)
  expect_equal(A, z)
})

test_that("blas2 dgerT", {
  alpha <- runif(1)
  x <- runif(2, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  z <- alpha * (matrix(y, ncol=1) %*% matrix(x, nrow=1)) + A
  C_dgerT(alpha, x, y, A)
  expect_equal(A, z)
})
