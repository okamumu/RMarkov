library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for GSSTEP")

# this test requires src/test_gsstep.cpp

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

test_that("gsstep vecN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- runif(4)
  x0 <- runif(4)
  res <- C_gsstep_vecN_dense(alpha, Q, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.vector(solve(D/omega + L) %*% (b/alpha - as.vector((U - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

test_that("gsstep vecT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- runif(4)
  x0 <- runif(4)
  res <- C_gsstep_vecT_dense(alpha, Q, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.vector(solve(D/omega + t(U)) %*% (b/alpha - as.vector((t(L) - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

test_that("gsstep vecN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- runif(4)
  x0 <- runif(4)
  res <- C_gsstep_vecN_csr(alpha, spQ, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.vector(solve(D/omega + L) %*% (b/alpha - as.vector((U - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

test_that("gsstep vecT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- runif(4)
  x0 <- runif(4)
  res <- C_gsstep_vecT_csc(alpha, spQ, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.vector(solve(D/omega + t(U)) %*% (b/alpha - as.vector((t(L) - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

#/// mat

test_that("gsstep matN dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- matrix(runif(4*3), 4, 3)
  x0 <- matrix(runif(4*3), 4, 3)
  res <- C_gsstep_matN_dense(alpha, Q, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.matrix(solve(D/omega + L) %*% (b/alpha - ((U - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

test_that("gsstep matT dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- matrix(runif(4*3), 4, 3)
  x0 <- matrix(runif(4*3), 4, 3)
  res <- C_gsstep_matT_dense(alpha, Q, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.matrix(solve(D/omega + t(U)) %*% (b/alpha - (t(L) - D * (1-omega)/omega - sigma * diag(4)) %*% x0))

  expect_equal(res, v)
})

test_that("gsstep matN csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- matrix(runif(4*3), 4, 3)
  x0 <- matrix(runif(4*3), 4, 3)
  res <- C_gsstep_matN_csr(alpha, spQ, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.matrix(solve(D/omega + L) %*% (b/alpha - ((U - D * (1-omega)/omega - sigma * diag(4)) %*% x0)))

  expect_equal(res, v)
})

test_that("gsstep matT csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  alpha <- runif(1)
  omega <- runif(1)
  sigma <- runif(1)
  b <- matrix(runif(4*3), 4, 3)
  x0 <- matrix(runif(4*3), 4, 3)
  res <- C_gsstep_matT_csc(alpha, spQ, sigma, omega, b, x0)

  D <- diag(diag(Q))
  U <- Q
  U[lower.tri(U, diag=T)] <- 0
  L <- Q
  L[upper.tri(L, diag=T)] <- 0
  v <- as.matrix(solve(D/omega + t(U)) %*% (b/alpha - (t(L) - D * (1-omega)/omega - sigma * diag(4)) %*% x0))

  expect_equal(res, v)
})
