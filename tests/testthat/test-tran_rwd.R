library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for tran rwd")

cumexp <- function(Q, t) {
  n <- dim(Q)[1]
  im <- diag(n)
  zm <- matrix(0,n,n)
  Qdash <- rbind(cbind(Q, im), cbind(zm, zm))
  expm(Qdash*t)[1:n,(n+1):(n+n)]
}

test_that("tran rwd vecT vecT forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(TRUE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(FALSE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT forward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT backward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT forward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT backward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT forward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT vecT backward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- c(1,2,3,4)
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_vec(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

########

test_that("tran rwd vecT matT forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(TRUE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(FALSE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT forward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT backward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT forward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT backward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT forward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(Q*sum(t)))
  cx <- as.vector(x0 %*% cumexp(Q, sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(Q*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(Q, t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd vecT matT backward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  cx0 <- c(0,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_vec_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(2, length(t)))
  res$crwd <- array(res$crwd, dim = c(2, length(t)))

  ct <- cumsum(t)
  x <- as.vector(x0 %*% expm(t(Q)*sum(t)))
  cx <- as.vector(x0 %*% cumexp(t(Q), sum(t)))
  irwd <- sapply(ct, function(t) as.vector(x0 %*% expm(t(Q)*t)) %*% rwd)
  crwd <- sapply(ct, function(t) as.vector(x0 %*% cumexp(t(Q), t)) %*% rwd)

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

########

test_that("tran rwd matT matT forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(TRUE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(t(Q)*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(t(Q), sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(t(Q)*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(t(Q), t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(FALSE, Matrix(Q), x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(Q*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(Q, sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(Q*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(Q, t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT forward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(t(Q)*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(t(Q), sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(t(Q)*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(t(Q), t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT backward csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(Q*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(Q, sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(Q*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(Q, t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT forward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(t(Q)*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(t(Q), sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(t(Q)*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(t(Q), t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT backward csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(Q*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(Q, sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(Q*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(Q, t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT forward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(TRUE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(t(Q)*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(t(Q), sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(t(Q)*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(t(Q), t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("tran rwd matT matT backward coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  cx0 <- matrix(0, 4, 2)
  rwd <- cbind(c(1,2,3,4), c(3,4,5,6), c(0,1,0,0))
  t <- c(0.3,0.4,0.6,0.1,0.2)
  res <- Ctran_unif_rwd_mat_mat(FALSE, spQ, x0, cx0, rwd, t, 1.01, 1.0e-8, 100)
  res$irwd <- array(res$irwd, dim = c(3, 2, length(t)))
  res$crwd <- array(res$crwd, dim = c(3, 2, length(t)))

  ct <- cumsum(t)
  x <- matrix(as.vector(expm(Q*sum(t)) %*% x0), 4, 2)
  cx <- matrix(as.vector(cumexp(Q, sum(t)) %*% x0), 4, 2)
  irwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% expm(Q*t) %*% x0))), dim=c(3,2,length(t)))
  crwd <- array(as.vector(sapply(ct, function(t) as.vector(t(rwd) %*% cumexp(Q, t) %*% x0))), dim=c(3,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

