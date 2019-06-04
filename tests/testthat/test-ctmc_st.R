library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for stationary analysis of CTMC")

# this test requires src/test_ctmc_st.cpp

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

test_that("gth dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- C_gth_dense(Q)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res, result)
})

test_that("gth csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "RsparseMatrix")
  res <- C_gth_csr(spQ)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res, result)
})

test_that("gth csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "CsparseMatrix")
  res <- C_gth_csc(spQ)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res, result)
})

test_that("gth coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  spQ <- as(Matrix(Q), "TsparseMatrix")
  res <- C_gth_coo(spQ)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res, result)
})

###

test_that("ctmc st gth dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.gth(Q)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res, result)
})

test_that("ctmc st gs dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.gs(Q=Q, matrix.class = "CsparseMatrix")
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res$x, result)
})

test_that("ctmc st power dense", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.power(Q=Q)
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res$x, result)
})

test_that("ctmc st power csr", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.power(Q=Q, matrix.class = "RsparseMatrix")
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res$x, result)
})

test_that("ctmc st power csc", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.power(Q=Q, matrix.class = "CsparseMatrix")
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res$x, result)
})

test_that("ctmc st power coo", {
  Q <- rbind(
    c(-9,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.st.power(Q=Q, matrix.class = "TsparseMatrix")
  A <- Q
  A[,1] <- 1
  result <- solve(t(A), c(1,0,0,0))
  expect_equal(res$x, result)
})
