library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for sensitivity analysis of stationary analysis of CTMC")

test_that("stsen gs csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  dQ <- rbind(
    c(-1,1,0,0),
    c(0,0,0,0),
    c(0,0,0,0),
    c(0,0,0,0)
  )
  pis <- ctmc.st.gth(Q)
  y <- as.vector(pis %*% dQ)
  x0 <- rep(1/4, 4)
  spQ <- as(Q, "CsparseMatrix")
  res <- Cmarkovstsen_gs(spQ, x0, y, pis, 50, 1.0e-7, 1000)

  A <- Q
  A[,1] <- 1
  y[1] <- 0
  result <- solve(t(A), -y)
  expect_equal(res$x, result)
})

test_that("stsen gs dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  dQ <- rbind(
    c(-1,1,0,0),
    c(0,0,0,0),
    c(0,0,0,0),
    c(0,0,0,0)
  )
  pis <- ctmc.st.gth(Q)
  y <- as.vector(pis %*% dQ)
  x0 <- rep(1/4, 4)
  spQ <- Matrix(Q)
  res <- Cmarkovstsen_gs(spQ, x0, y, pis, 50, 1.0e-7, 1000)

  A <- Q
  A[,1] <- 1
  y[1] <- 0
  result <- solve(t(A), -y)
  expect_equal(res$x, result)
})

##

test_that("ctmc.stsen.gs", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  dQ <- rbind(
    c(-1,1,0,0),
    c(0,0,0,0),
    c(0,0,0,0),
    c(0,0,0,0)
  )
  pis <- ctmc.st.gth(Q)
  res <- ctmc.stsen.gs(Q=Q, pis=pis, b=pis %*% dQ)

  A <- Q
  A[,1] <- 1
  y <- as.vector(pis %*% dQ)
  y[1] <- 0
  result <- solve(t(A), -y)
  expect_equal(res$x, result)
})

