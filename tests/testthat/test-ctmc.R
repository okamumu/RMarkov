library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for utilities of CTMC")

# this test requires src/test_ctmc.cpp

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

test_that("unif dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-20,7,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res = C_unif_dense(Q, 1.01)
  P <- diag(4) + Q / (max(abs(diag(Q))) * 1.01)
  expect_equal(res$P, P)
})

test_that("unif csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-20,7,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  res = C_unif_csr(spQ, 1.01)
  P <- as(Matrix(diag(4) + Q / (max(abs(diag(Q))) * 1.01)), "RsparseMatrix")
  expect_equal(res$P@x, P@x)
  expect_equal(res$P@p, P@p)
  expect_equal(res$P@j, P@j)
})

test_that("unif csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-20,7,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  res = C_unif_csc(spQ, 1.01)
  P <- as(Matrix(diag(4) + Q / (max(abs(diag(Q))) * 1.01)), "CsparseMatrix")
  expect_equal(res$P@x, P@x)
  expect_equal(res$P@p, P@p)
  expect_equal(res$P@i, P@i)
})

test_that("unif coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-20,7,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  res = C_unif_coo(spQ, 1.01)
  P <- as(Matrix(diag(4) + Q / (max(abs(diag(Q))) * 1.01)), "TsparseMatrix")
  expect_equal(res$P@x, P@x)
  expect_equal(res$P@i, P@i)
  expect_equal(res$P@j, P@j)
})

## mexp

test_that("mexp vecT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  y <- C_mexp_vecT_dense(Q, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(x %*% expm(Q*t)))
})

test_that("mexp vecN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  y <- C_mexp_vecN_dense(Q, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(expm(Q*t) %*% x))
})

test_that("mexp vecT csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  y <- C_mexp_vecT_csr(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(x %*% expm(Q*t)))
})

test_that("mexp vecN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  y <- C_mexp_vecN_csr(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(expm(Q*t) %*% x))
})

test_that("mexp vecT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  y <- C_mexp_vecT_csc(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(x %*% expm(Q*t)))
})

test_that("mexp vecN csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  y <- C_mexp_vecN_csc(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(expm(Q*t) %*% x))
})

test_that("mexp vecT coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  y <- C_mexp_vecT_coo(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(x %*% expm(Q*t)))
})

test_that("mexp vecN coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  y <- C_mexp_vecN_coo(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expect_equal(y, as.vector(expm(Q*t) %*% x))
})

## mexp

test_that("mexp matN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matN_dense(Q, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(Q*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matT_dense(Q, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(t(Q)*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matN_csr(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(Q*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matT csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matT_csr(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(t(Q)*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matN csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matN_csc(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(Q*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matT_csc(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(t(Q)*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matN coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matN_coo(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(Q*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

test_that("mexp matT coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4*3), 4, 3)
  y <- C_mexp_matT_coo(spQ, t, x, eps = 1.0e-8, ufactor = 1.01)
  expected <- as.matrix(expm(t(Q)*t) %*% x)
  attr(expected, "dimnames") <- NULL
  expect_equal(y, expected)
})

## mexpint

test_that("mexpint vecT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  cy <- runif(4)
  res <- C_mexpint_vecT_dense(Q, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(x %*% expm(Q*t))
  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  ex.cy <- as.vector(c(x, rep(0,4)) %*% expm(Qdash*t))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  cy <- runif(4)
  res <- C_mexpint_vecN_dense(Q, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(expm(Q*t) %*% x)
  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  ex.cy <- as.vector(expm(Qdash*t) %*% c(x, rep(0,4)))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecT csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  cy <- runif(4)
  res <- C_mexpint_vecT_csr(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(x %*% expm(Q*t))
  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  ex.cy <- as.vector(c(x, rep(0,4)) %*% expm(Qdash*t))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  cy <- runif(4)
  res <- C_mexpint_vecN_csr(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(expm(Q*t) %*% x)
  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  ex.cy <- as.vector(expm(Qdash*t) %*% c(x, rep(0,4)))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  cy <- runif(4)
  res <- C_mexpint_vecT_csc(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(x %*% expm(Q*t))
  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  ex.cy <- as.vector(c(x, rep(0,4)) %*% expm(Qdash*t))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecN csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  cy <- runif(4)
  res <- C_mexpint_vecN_csc(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(expm(Q*t) %*% x)
  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  ex.cy <- as.vector(expm(Qdash*t) %*% c(x, rep(0,4)))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecT coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  x <- x / sum(x)
  cy <- runif(4)
  res <- C_mexpint_vecT_coo(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(x %*% expm(Q*t))
  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  ex.cy <- as.vector(c(x, rep(0,4)) %*% expm(Qdash*t))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

test_that("mexpint vecN coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- runif(4)
  cy <- runif(4)
  res <- C_mexpint_vecN_coo(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)
  ex.y <- as.vector(expm(Q*t) %*% x)
  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  ex.cy <- as.vector(expm(Qdash*t) %*% c(x, rep(0,4)))[5:8] + cy
  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy)
})

## mexpint mat

test_that("mexpint matT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matT_dense(Q, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(t(Q)*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(t(Qdash)*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matN_dense(Q, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(Q*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(Qdash*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matT csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matT_csr(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(t(Q)*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(t(Qdash)*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matN_csr(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(Q*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(Qdash*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matT_csc(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(t(Q)*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(t(Qdash)*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matN csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matN_csc(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(Q*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(Qdash*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matT coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matT_coo(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(t(Q)*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- rbind(cbind(Q, diag(4)), matrix(0, 4, 8))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(t(Qdash)*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

test_that("mexpint matN coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  t <- runif(1, min=0.1, max=1.0)
  x <- matrix(runif(4), 4, 3)
  cy <- matrix(runif(4), 4, 3)
  cy <- matrix(0, 4, 3)
  res <- C_mexpint_matN_coo(spQ, t, x, cy, eps = 1.0e-8, ufactor = 1.01)

  ex.y <- as.matrix(expm(Q*t) %*% x)
  attr(ex.y, "dimnames") <- NULL

  Qdash <- cbind(rbind(Q, diag(4)), matrix(0, 8, 4))
  xdash <- rbind(x, matrix(0,4,3))
  ex.cy <- as.matrix(expm(Qdash*t) %*% xdash)[5:8,]
  attr(ex.cy, "dimnames") <- NULL

  expect_equal(res$y, ex.y)
  expect_equal(res$cy, ex.cy + cy)
})

