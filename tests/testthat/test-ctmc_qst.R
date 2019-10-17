library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for quasi-stationary analysis of CTMC")

# this test requires src/test_ctmc_qst.cpp

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

test_that("ctmc st gs csc", {
  Q <- rbind(
    c(-15,2,3,4),
    c(5,-13,0,8),
    c(9,1,-13,3),
    c(5,2,1,-8)
  )
  res <- ctmc.qst.gs(Q=Q, matrix.class = "CsparseMatrix")
  result <- ctmc.st.power(Q=Q, matrix.class = "CsparseMatrix")
  expect_equal(res$x, result$x)
})
