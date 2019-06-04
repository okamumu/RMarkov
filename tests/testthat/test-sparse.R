library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for utilities of sparse")

# this test requires src/test_sparse.cpp

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

test_that("sparse nnz", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  expect_equal(C_sparse_nnz(A), sum(A!=0))
})

test_that("nrow, ncol dense", {
  m <- sample(x=1:100, size=1)
  n <- sample(x=1:100, size=1)
  A <- matrix(runif(m*n), m, n)
  d <- C_nrow_dense(A)
  expect_equal(d, c(m,n))
})

test_that("nrow, ncol csr", {
  m <- sample(x=1:100, size=1)
  n <- sample(x=1:100, size=1)
  A <- matrix(runif(m*n), m, n)
  spA <- as(Matrix(A), "RsparseMatrix")
  d <- C_nrow_csr(spA)
  expect_equal(d, c(m,n))
})

test_that("nrow, ncol csc", {
  m <- sample(x=1:100, size=1)
  n <- sample(x=1:100, size=1)
  A <- matrix(runif(m*n), m, n)
  spA <- as(Matrix(A), "CsparseMatrix")
  d <- C_nrow_csc(spA)
  expect_equal(d, c(m,n))
})

test_that("nrow, ncol coo", {
  m <- sample(x=1:100, size=1)
  n <- sample(x=1:100, size=1)
  A <- matrix(runif(m*n), m, n)
  spA <- as(Matrix(A), "TsparseMatrix")
  d <- C_nrow_coo(spA)
  expect_equal(d, c(m,n))
})

test_that("dense to sparse CSR", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  nnz <- C_sparse_nnz(A)
  value <- numeric(nnz)
  ptr <- integer(nrow(A)+1)
  ind <- integer(nnz)
  C_create_csr(A, value, ptr, ind)
  B <- as(Matrix(A), "RsparseMatrix")
  expect_equal(value, B@x)
  expect_equal(ptr, B@p)
  expect_equal(ind, B@j)
})

test_that("sparse to dense CSR", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  A <- as(Matrix(A), "RsparseMatrix")
  B <- matrix(0, 10, 20)
  C <- as.matrix(A)
  C_sparse_to_dense_csr(A, B)
  attr(C, "dimnames") <- NULL
  expect_equal(B, C)
})

test_that("dense to sparse CSC", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  nnz <- C_sparse_nnz(A)
  value <- numeric(nnz)
  ptr <- integer(ncol(A)+1)
  ind <- integer(nnz)
  C_create_csc(A, value, ptr, ind)
  B <- as(Matrix(A), "CsparseMatrix")
  expect_equal(value, B@x)
  expect_equal(ptr, B@p)
  expect_equal(ind, B@i)
})

test_that("sparse to dense CSC", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  A <- as(Matrix(A), "CsparseMatrix")
  B <- matrix(0, 10, 20)
  C <- as.matrix(A)
  C_sparse_to_dense_csc(A, B)
  attr(C, "dimnames") <- NULL
  expect_equal(B, C)
})

test_that("dense to sparse COO", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  nnz <- C_sparse_nnz(A)
  value <- numeric(nnz)
  row <- integer(nnz)
  col <- integer(nnz)
  C_create_coo(A, value, row, col)
  B <- as(Matrix(A), "TsparseMatrix")
  expect_equal(value, B@x)
  expect_equal(row, B@i)
  expect_equal(col, B@j)
})

test_that("sparse to dense COO", {
  p <- runif(1, min=0.1, max=0.9)
  A <- create.sparse(10, 20, p)
  A <- as(Matrix(A), "TsparseMatrix")
  B <- matrix(0, 10, 20)
  C <- as.matrix(A)
  C_sparse_to_dense_coo(A, B)
  attr(C, "dimnames") <- NULL
  expect_equal(B, C)
})

test_that("diag_get dense", {
  m <- 10
  n <- 5
  A <- matrix(runif(m*n), m, n)
  x <- C_diag_dense(A)
  expect_equal(x, diag(A))
})

test_that("diag_get csr", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
    )
  spA <- as(Matrix(A), "RsparseMatrix")
  x <- C_diag_csr(spA)
  expect_equal(x, diag(A))
})

test_that("diag_get csr (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "RsparseMatrix")
  x <- C_diag_csr(spA)
  expect_equal(x, diag(A))
})

test_that("diag_get csc", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "CsparseMatrix")
  x <- C_diag_csc(spA)
  expect_equal(x, diag(A))
})

test_that("diag_get csc (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "CsparseMatrix")
  x <- C_diag_csc(spA)
  expect_equal(x, diag(A))
})

test_that("diag_get coo", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "TsparseMatrix")
  x <- C_diag_coo(spA)
  expect_equal(x, diag(A))
})

test_that("diag_get coo (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "TsparseMatrix")
  x <- C_diag_coo(spA)
  expect_equal(x, diag(A))
})


### diag set

test_that("diag_set dense", {
  m <- 10
  n <- 5
  A <- matrix(runif(m*n), m, n)
  x <- runif(min(m,n))
  C_diag_set_dense(A, x)
  expect_equal(diag(A), x)
})

test_that("diag_set csr", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "RsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_csr(spA, x)
  expect_equal(diag(as.matrix(spA)), x)
})

test_that("diag_set csr (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  b <- diag(A) != 0
  spA <- as(Matrix(A), "RsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_csr(spA, x)
  expect_equal(all(diag(as.matrix(spA))[!b]==0), TRUE)
  expect_equal(diag(as.matrix(spA))[b], x[b])
})

test_that("diag_set csc", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "CsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_csc(spA, x)
  expect_equal(diag(as.matrix(spA)), x)
})

test_that("diag_set csc (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  b <- diag(A) != 0
  spA <- as(Matrix(A), "CsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_csc(spA, x)
  expect_equal(all(diag(as.matrix(spA))[!b]==0), TRUE)
  expect_equal(diag(as.matrix(spA))[b], x[b])
})

test_that("diag_set coo", {
  A <- rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,1,2,3)
  )
  spA <- as(Matrix(A), "TsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_coo(spA, x)
  expect_equal(diag(as.matrix(spA)), x)
})

test_that("diag_set coo (diag elements include 0)", {
  A <- rbind(
    c(1,2,3,4),
    c(5,0,7,8),
    c(9,1,2,3)
  )
  b <- diag(A) != 0
  spA <- as(Matrix(A), "TsparseMatrix")
  x <- runif(min(nrow(A), ncol(A)))
  C_diag_set_coo(spA, x)
  expect_equal(all(diag(as.matrix(spA))[!b]==0), TRUE)
  expect_equal(diag(as.matrix(spA))[b], x[b])
})
