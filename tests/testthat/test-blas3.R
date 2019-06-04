library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for blas3 (dense matrix)")

# this test requires src/test_blas3.cpp

test_that("blas3 dgemmNN", {
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(2, min = -10, max = 10), 2, 4)
  C <- matrix(runif(3, min = -10, max = 10), 3, 4)
  z <- alpha * (A %*% B) + beta * C
  C_dgemmNN(alpha, A, B, beta, C)
  expect_equal(C, z)
})

test_that("blas3 dgemmTN", {
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(6, min = -10, max = 10), 2, 3)
  B <- matrix(runif(2, min = -10, max = 10), 2, 4)
  C <- matrix(runif(3, min = -10, max = 10), 3, 4)
  z <- alpha * (t(A) %*% B) + beta * C
  C_dgemmTN(alpha, A, B, beta, C)
  expect_equal(C, z)
})

test_that("blas3 dgemmNT", {
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(2, min = -10, max = 10), 4, 2)
  C <- matrix(runif(3, min = -10, max = 10), 3, 4)
  z <- alpha * (A %*% t(B)) + beta * C
  C_dgemmNT(alpha, A, B, beta, C)
  expect_equal(C, z)
})

test_that("blas3 dgemmTT", {
  alpha <- runif(1)
  beta <- runif(1)
  A <- matrix(runif(6, min = -10, max = 10), 2, 3)
  B <- matrix(runif(2, min = -10, max = 10), 4, 2)
  C <- matrix(runif(3, min = -10, max = 10), 3, 4)
  z <- alpha * (t(A) %*% t(B)) + beta * C
  C_dgemmTT(alpha, A, B, beta, C)
  expect_equal(C, z)
})
