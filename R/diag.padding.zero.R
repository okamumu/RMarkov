#' Padding diagonal elements
#'
#' The function executes the padding diagonal elements as zero if they are skipped in the sparse matrix representation.
#'
#' @param A A  matrix. This is one of matrix, dgeMatrix, dgCMatrix, dgRMatrix and dgTMatrix

diag.padding.zero <- function(A) {
  cls <- class(A)
  x <- Matrix::diag(A)
  Matrix::diag(A)[x == 0] <- NA
  A@x[is.na(A@x)] <- 0
  switch(cls,
         "dgeMatrix" = as(A, "dgeMatrix"),
         "dgRMatrix" = as(A, "RsparseMatrix"),
         "dgCMatrix" = as(A, "CsparseMatrix"),
         "dgTMatrix" = as(A, "TsparseMatrix"),
         stop("invalid Matrix class:", cls)
  )
}
