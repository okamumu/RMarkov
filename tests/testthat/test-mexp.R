library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for mexp")

test_that("mexp vec", {
  A <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  y <- mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100)
  expect_equal(y, as.vector(expm(A*t) %*% x))
  y <- mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100)
  expect_equal(y, as.vector(expm(t(A)*t) %*% x))
  y <- mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "RsparseMatrix")
  expect_equal(y, as.vector(expm(A*t) %*% x))
  y <- mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "RsparseMatrix")
  expect_equal(y, as.vector(expm(t(A)*t) %*% x))
  y <- mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "CsparseMatrix")
  expect_equal(y, as.vector(expm(A*t) %*% x))
  y <- mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "CsparseMatrix")
  expect_equal(y, as.vector(expm(t(A)*t) %*% x))
  y <- mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "TsparseMatrix")
  expect_equal(y, as.vector(expm(A*t) %*% x))
  y <- mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "TsparseMatrix")
  expect_equal(y, as.vector(expm(t(A)*t) %*% x))
})

test_that("mexp mat", {
  A <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  x <- x / sum(x)
  y <- as(mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100), "dgeMatrix")
  expect_equal(y, expm(A*t) %*% x)
  y <- as(mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100), "dgeMatrix")
  expect_equal(y, expm(t(A)*t) %*% x)
  y <- as(mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "RsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(A*t) %*% x)
  y <- as(mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "RsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(t(A)*t) %*% x)
  y <- as(mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "CsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(A*t) %*% x)
  y <- as(mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "CsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(t(A)*t) %*% x)
  y <- as(mexpAx(A, x, t, FALSE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "TsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(A*t) %*% x)
  y <- as(mexpAx(A, x, t, TRUE, eps = 1.0e-8, ufact = 1.01, rmax = 100, matrix.class = "TsparseMatrix"), "dgeMatrix")
  expect_equal(y, expm(t(A)*t) %*% x)
})

