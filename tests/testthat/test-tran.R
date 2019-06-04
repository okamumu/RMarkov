library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for tran")

cumexp <- function(Q, t) {
  n <- dim(Q)[1]
  im <- diag(n)
  zm <- matrix(0,n,n)
  Qdash <- rbind(cbind(Q, im), cbind(zm, zm))
  expm(Qdash*t)[1:n,(n+1):(n+n)]
}

test_that("mexp vecT forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(TRUE, Matrix(Q), x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)))
  cx <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(FALSE, Matrix(Q), x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})


test_that("mexp vecT forward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)))
  cx <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT backward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT forward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)))
  cx <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT backward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT forward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)))
  cx <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp vecT backward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_vec(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,length(t)))
  res$cx <- array(res$cx, dim=c(4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

#####

test_that("mexp matT forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(TRUE, Matrix(Q), x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(t(Q)*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(t(Q), t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(FALSE, Matrix(Q), x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT forward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(t(Q)*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(t(Q), t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT backward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT forward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(t(Q)*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(t(Q), t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT backward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT forward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(TRUE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(t(Q)*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(t(Q), t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("mexp matT backward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- rbind(
    c(1,1,0,0),
    c(0,1,0,0),
    c(1,0,1,0),
    c(0,0,0,2)
  )
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_mat(FALSE, spQ, x0, t, 1.01, 1.0e-8, 100)
  res$x <- array(res$x, dim=c(4,4,length(t)))
  res$cx <- array(res$cx, dim=c(4,4,length(t)))

  ct <- cumsum(t)
  x <- sapply(ct, function(t) as.vector(expm(Q*t) %*% x0))
  x <- array(x, dim=c(4,4,length(t)))
  cx <- sapply(ct, function(t) as.vector(cumexp(Q, t) %*% x0))
  cx <- array(cx, dim=c(4,4,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})
