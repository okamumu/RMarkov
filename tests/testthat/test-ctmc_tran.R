library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for ctmc.tran")

cumexp <- function(Q, t) {
  n <- dim(Q)[1]
  im <- diag(n)
  zm <- matrix(0,n,n)
  Qdash <- rbind(cbind(Q, im), cbind(zm, zm))
  expm(Qdash*t)[1:n,(n+1):(n+n)]
}

test_that("ctmc tran vec forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.unif(Q, x0, t)

  x <- sapply(t, function(u) as.vector(expm(t(Q)*u) %*% x0))
  cx <- sapply(t, function(u) as.vector(cumexp(t(Q), u) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("ctmc tran vec backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.unif(Q, x0, t, transpose = FALSE)

  x <- sapply(t, function(u) as.vector(expm(Q*u) %*% x0))
  cx <- sapply(t, function(u) as.vector(cumexp(Q, u) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("ctmc tran vec forward sparse", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.unif(Q, x0, t, matrix.class = "RsparseMatrix")

  x <- sapply(t, function(u) as.vector(expm(t(Q)*u) %*% x0))
  cx <- sapply(t, function(u) as.vector(cumexp(t(Q), u) %*% x0))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

test_that("ctmc tran mat forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0), c(0,0,1,0))
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.unif(Q, x0, t, matrix.class = "RsparseMatrix")

  x <- as.vector(sapply(t, function(u) as.vector(expm(t(Q)*u) %*% x0)))
  x <- array(x, dim = c(4,2,length(t)))
  cx <- as.vector(sapply(t, function(u) as.vector(cumexp(t(Q), u) %*% x0)))
  cx <- array(cx, dim = c(4,2,length(t)))

  expect_equal(res$x, x)
  expect_equal(res$cx, cx)
})

###

test_that("ctmc tran rwd vec vec forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  rwd <- c(1,2,3,4)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t)

  irwd <- sapply(t, function(u) as.vector(rwd %*% expm(t(Q)*u) %*% x0))
  crwd <- sapply(t, function(u) as.vector(rwd %*% cumexp(t(Q), u) %*% x0))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("ctmc tran rwd vec vec backward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  rwd <- c(1,2,3,4)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t, transpose = FALSE)

  irwd <- sapply(t, function(u) as.vector(rwd %*% expm(Q*u) %*% x0))
  crwd <- sapply(t, function(u) as.vector(rwd %*% cumexp(Q, u) %*% x0))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("ctmc tran rwd vec vec forward sparse", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  rwd <- c(1,2,3,4)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t, matrix.class = "RsparseMatrix")

  irwd <- sapply(t, function(u) as.vector(rwd %*% expm(t(Q)*u) %*% x0))
  crwd <- sapply(t, function(u) as.vector(rwd %*% cumexp(t(Q), u) %*% x0))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("ctmc tran rwd vec mat forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- c(1,0,0,0)
  rwd <- cbind(c(1,2,3,4), c(5,6,7,8))
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t)
  irwd <- sapply(t, function(u) as.vector(t(rwd) %*% expm(t(Q)*u) %*% x0))
  crwd <- sapply(t, function(u) as.vector(t(rwd) %*% cumexp(t(Q), u) %*% x0))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("ctmc tran rwd mat mat forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  rwd <- cbind(c(1,2,3,4), c(5,6,7,8), c(4,5,6,7))
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t)

  irwd <- as.vector(sapply(t, function(u) as.vector(t(rwd) %*% expm(t(Q)*u) %*% x0)))
  crwd <- as.vector(sapply(t, function(u) as.vector(t(rwd) %*% cumexp(t(Q), u) %*% x0)))
  irwd <- array(irwd, dim = c(3,2,length(t)))
  crwd <- array(crwd, dim = c(3,2,length(t)))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})

test_that("ctmc tran rwd mat vec forward dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,0,1,-6)
  )
  x0 <- cbind(c(1,0,0,0),c(0,1,0,0))
  rwd <- c(1,2,3,4)
  t <- seq(1,10,length.out = 100)
  res <- ctmc.tran.rwd.unif(Q, x0, rwd, t)

  irwd <- as.vector(sapply(t, function(u) as.vector(rwd %*% expm(t(Q)*u) %*% x0)))
  crwd <- as.vector(sapply(t, function(u) as.vector(rwd %*% cumexp(t(Q), u) %*% x0)))
  irwd <- array(irwd, dim = c(1,2,length(t)))
  crwd <- array(crwd, dim = c(1,2,length(t)))

  expect_equal(res$irwd, irwd)
  expect_equal(res$crwd, crwd)
})
