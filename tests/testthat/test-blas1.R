library(testthat)
library(RMarkov)
library(Matrix)
library(Rcpp)

context("Unit tests for blas1")

x <- runif(3, min = -10, max = 10)
y <- runif(3, min = -10, max = 10)

A <- matrix(runif(6, min = -10, max = 10), 3, 2)
B <- matrix(runif(6, min = -10, max = 10), 3, 2)
C <- matrix(runif(6, min = -10, max = 10), 3, 2)

geA <- Matrix(A)
geB <- Matrix(B)

test_that("blas1 ddot", {
  expect_equal(C_ddot(x,y), sum(x*y))
  expect_equal(C_ddot(A,B), sum(A*B))
  expect_equal(C_ddot2(geA,geB), sum(A*B))
})

test_that("blas1 dasum", {
  expect_equal(C_dasum(x), sum(abs(x)))
  expect_equal(C_dasum(A), sum(abs(A)))
  expect_equal(C_dasum2(geA), sum(abs(A)))
})

test_that("blas1 idamax", {
  expect_equal(C_idamax(x), which.max(abs(x))-1)
  expect_equal(C_idamax(A), which.max(abs(A))-1)
  expect_equal(C_idamax2(geA), which.max(abs(A))-1)
})

test_that("blas1 dscal", {
  alpha <- runif(1)
  C_dcopy(x, y)
  C_dscal(alpha, x)
  expect_equal(x, y * alpha)
  C_dcopy(A, B)
  C_dscal(alpha, A)
  expect_equal(A, B * alpha)
})

test_that("blas1 dcopy", {
  C_dcopy(x, y)
  expect_equal(y, x)
  C_dcopy(A, B)
  expect_equal(B, A)
  C_dcopy2(geA, C)
  expect_equal(as.vector(C), as.vector(as.matrix(geA)))
})

