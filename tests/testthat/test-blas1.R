library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for blas1")

# this test requires src/test_blas1.cpp

test_that("blas1 ddot", {
  x <- runif(3, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(6, min = -10, max = 10), 3, 2)
  geA <- Matrix(A)
  expect_equal(C_ddot(x,y), sum(x*y))
  expect_equal(C_ddot(A,B), sum(A*B))
  expect_equal(C_ddot2(geA,B), sum(A*B)) # S4
})

test_that("blas1 dasum", {
  x <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  expect_equal(C_dasum(x), sum(abs(x)))
  expect_equal(C_dasum(A), sum(abs(A)))
})

test_that("blas1 idamax", {
  x <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  expect_equal(C_idamax(x), which.max(abs(x))-1)
  expect_equal(C_idamax(A), which.max(abs(A))-1)
})

test_that("blas1 dscal", {
  alpha <- runif(1)
  x <- runif(3, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(6, min = -10, max = 10), 3, 2)
  C_dcopy(x, y)
  C_dscal(alpha, x)
  expect_equal(x, y * alpha)
  C_dcopy(A, B)
  C_dscal(alpha, A)
  expect_equal(A, B * alpha)
})

test_that("blas1 daxpy", {
  alpha <- runif(1)
  x <- runif(3, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(6, min = -10, max = 10), 3, 2)
  z <- alpha * x + y
  C_daxpy(alpha, x, y)
  expect_equal(y, z)
  C <- alpha * A + B
  C_daxpy(alpha, A, B)
  expect_equal(B, C)
})

test_that("blas1 dcopy", {
  x <- runif(3, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(6, min = -10, max = 10), 3, 2)
  geA <- Matrix(A)
  C_dcopy(x, y)
  expect_equal(y, x)
  C_dcopy(A, B)
  expect_equal(B, A)
  C_dcopy2(geA, B)
  expect_equal(as.vector(B), as.vector(as.matrix(geA)))
})

test_that("blas1 dfill", {
  v <- runif(1)
  x <- runif(3, min = -10, max = 10)
  y <- runif(3, min = -10, max = 10)
  A <- matrix(runif(6, min = -10, max = 10), 3, 2)
  B <- matrix(runif(6, min = -10, max = 10), 3, 2)
  C_dfill(x, v)
  expect_equal(x, rep(v, 3))
  C_dfill(A, v)
  expect_equal(A, matrix(rep(v, 6), 3, 2))
})

