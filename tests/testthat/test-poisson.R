library(testthat)
library(RMarkov)
library(Matrix)

context("Unit tests for poisson")

# this test requires src/test_poisson.cpp

test_that("poisson rightbound 1.0e-6", {
  eps <- 1.0e-6
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-7", {
  eps <- 1.0e-7
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-8", {
  eps <- 1.0e-8
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-9", {
  eps <- 1.0e-9
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-10", {
  eps <- 1.0e-10
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-11", {
  eps <- 1.0e-11
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson rightbound 1.0e-12", {
  eps <- 1.0e-12
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  r_e <- qpois(eps, lambda, lower.tail = FALSE)
  expect_gte(r, r_e)
})

test_that("poisson prob", {
  eps <- 1.0e-12
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  res <- C_poi_pmf(lambda, 0, r)
  res0 <- dpois(x=0:r, lambda = lambda)
  expect_equal(res$prob/res$weight, res0)
})

test_that("poisson cprob", {
  eps <- 1.0e-12
  lambda <- rexp(1, rate=1/1.0e+5)
  r <- C_poi_rightbound(lambda, eps)
  res <- C_poi_cpmf(lambda, 0, r)
  res0 <- ppois(q=0:r, lambda = lambda, lower.tail = FALSE)
  expect_equal(res$cprob/res$weight, res0)
})
