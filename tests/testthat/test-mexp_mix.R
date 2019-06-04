library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for mexp.mix")

cumexp <- function(Q, t) {
  n <- dim(Q)[1]
  im <- diag(n)
  zm <- matrix(0,n,n)
  Qdash <- rbind(cbind(Q, im), cbind(zm, zm))
  expm(Qdash*t)[1:n,(n+1):(n+n)]
}

#########

test_that("mexp mix vec notrans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(FALSE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)

})

test_that("mexp mix vec trans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(TRUE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)
})

test_that("mexp mix vec notrans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)

})

test_that("mexp mix vec trans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)
})

test_that("mexp mix vec notrans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)

})

test_that("mexp mix vec trans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)
})

test_that("mexp mix vec notrans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)

})

test_that("mexp mix vec trans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)

  expect_equal(res, x)
})

#########

test_that("mexp mix mat notrans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(FALSE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat trans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(TRUE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat notrans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat trans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat notrans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat trans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat notrans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

test_that("mexp mix mat trans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0), c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexp_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)

  expect_equal(res, x)
})

#########

test_that("mexpint mix vec notrans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(FALSE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec trans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(TRUE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec notrans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec trans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec notrans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec trans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec notrans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix vec trans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- c(1,0,0,0)
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_vec(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- apply(x,1,sum)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- apply(cx,1,sum)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

####

test_that("mexpint mix mat notrans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(FALSE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat trans dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(TRUE, Matrix(Q), x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat notrans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat trans csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "RsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat notrans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat trans csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "CsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat notrans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(FALSE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(Q*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(Q, t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

test_that("mexpint mix mat trans coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  spQ <- as(Q, "TsparseMatrix")
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  w <- runif(10)
  t <- seq(0,1,length.out = 10)
  dt <- diff(c(0,t))
  res <- Cmexpint_mix_unif_mat(TRUE, spQ, x0, w, dt, 1.01, 1.0e-8, 100)

  x <- sapply(1:length(t), function(i) w[i] * as.vector(expm(t(Q)*t[i]) %*% x0))
  x <- matrix(apply(x,1,sum), 4, 2)
  cx <- sapply(1:length(t), function(i) w[i] * as.vector(cumexp(t(Q), t[i]) %*% x0))
  cx <- matrix(apply(cx,1,sum), 4, 2)

  expect_equal(res$y, x)
  expect_equal(res$cy, cx)
})

#####

test_that("mexp mix vec integrate", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  res <- mexpAx.mix(Q, x=x0, f=dgamma, shape=2, scale=10)

  ff <- function(u, i) {
    as.vector(expm(Q*u) %*% x0)[i] * dgamma(u, shape=2, scale=10)
  }
  res2 <- sapply(1:4, function(i) deformula::deformula.zeroinf(f=ff, i=i)$value)
  expect_equal(res, res2)
})
