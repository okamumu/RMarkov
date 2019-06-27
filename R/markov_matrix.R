#
# StructuredMatrix <- R6::R6Class(
#   private = list(
#     dim = NA
#   ),
#   public = list(
#     initialize = function(m) {}
#   )
# )
#
# BlockMatrix <- R6::R6Class(
#   inherit = StructuredMatrix,
#   private = list(
#     blockDim = NA,
#     x = NA,
#     i = NA,
#     j = NA
#   ),
#   public = list(
#     initialize = function(i, j, x, dims, si, sj) {
#       if (missing(dims)) {
#         if (missing(si)) {
#           mi <- max(i)
#         } else {
#           mi <- length(si)
#         }
#         if (missing(sj)) {
#           mj <- max(j)
#         } else {
#           mj <- length(sj)
#         }
#         dims <- c(mi, mj)
#       }
#       if (missing(si)) {
#         si <- integer(dims[1])
#       }
#       if (missing(sj)) {
#         sj <- integer(dims[2])
#       }
#       stopifnot(length(si) == dims[1], length(sj) == dims[2])
#
#       if (length(x) != 0) {
#         for (u in 1:length(x)) {
#           v <- i[u]
#           si[v] <- dim(x[[u]])[1]
#         }
#         for (u in 1:length(x)) {
#           v <- j[u]
#           sj[v] <- dim(x[[u]])[2]
#         }
#       }
#
#       # check si sj
#       #  stopifnot(all(si != 0) && all(sj != 0))
#
#       mdims <- c(sum(si), sum(sj))
#       # ii <- cumsum(c(1,si))[1:dims[1]]
#       # jj <- cumsum(c(1,sj))[1:dims[2]]
#     },
#     set_name = function(name) { private$name <- name },
#     get_name = function() { private$name }
#   )
# )
#
# me <- Person$new("hoxo_m")
# me$get_name()
#
# setClass("sSparseMatrix", representation(MDim="integer"))
#
# setMethod("dim", signature(x = "sSparseMatrix"), function(x) x@MDim)
# setMethod("as.matrix", signature(x = "sSparseMatrix"), function(x) as(x, "TsparseMatrix"))
#
# #### blockMatrix
#
# setClass("bSparseMatrix", representation(i="integer", j="integer", Dim="integer",
#                                          x="list", si="integer", sj="integer"), contains="sSparseMatrix")
#
# blockMatrix <- function(i=c(), j=c(), x=list(), dims, si, sj) {
#
#   stopifnot(length(i) == length(j) && length(i) == length(x))
#
#   if (missing(dims)) {
#     if (missing(si)) {
#       mi <- max(i)
#     } else {
#       mi <- length(si)
#     }
#     if (missing(sj)) {
#       mj <- max(j)
#     } else {
#       mj <- length(sj)
#     }
#     dims <- c(mi, mj)
#   }
#   if (missing(si)) {
#     si <- integer(dims[1])
#   }
#   if (missing(sj)) {
#     sj <- integer(dims[2])
#   }
#   stopifnot(length(si) == dims[1], length(sj) == dims[2])
#
#   if (length(x) != 0) {
#     for (u in 1:length(x)) {
#       v <- i[u]
#       si[v] <- dim(x[[u]])[1]
#     }
#     for (u in 1:length(x)) {
#       v <- j[u]
#       sj[v] <- dim(x[[u]])[2]
#     }
#   }
#
#   # check si sj
#   #  stopifnot(all(si != 0) && all(sj != 0))
#
#   mdims <- c(sum(si), sum(sj))
#   # ii <- cumsum(c(1,si))[1:dims[1]]
#   # jj <- cumsum(c(1,sj))[1:dims[2]]
#
#   new("bSparseMatrix", i=as.integer(i), j=as.integer(j), Dim=as.integer(dims),
#       x=x, si=as.integer(si), sj=as.integer(sj), MDim=as.integer(mdims))
# }
#
# convBtoM <- function(bm) {
#   if (length(bm@x) == 0) {
#     return(Matrix(0, bm@MDim[1], bm@MDim[2]))
#   }
#   ii <- cumsum(c(1,bm@si))[1:bm@Dim[1]]
#   jj <- cumsum(c(1,bm@sj))[1:bm@Dim[2]]
#   ix <- c()
#   jx <- c()
#   xx <- c()
#   for (u in 1:length(bm@x)) {
#     i <- ii[bm@i[u]]
#     j <- jj[bm@j[u]]
#     y <- as(bm@x[[u]], "TsparseMatrix")
#     ix <- c(ix, i + y@i)
#     jx <- c(jx, j + y@j)
#     xx <- c(xx, y@x)
#   }
#   sparseMatrix(i=ix, j=jx, x=xx, dims=bm@MDim)
# }
#
# setMethod("[", signature(x = "bSparseMatrix", i = "index", j = "index", drop = "ANY"),
#           function (x, i, j, ..., drop) {
#             stopifnot(length(i) == 1 && length(j) == 1)
#             if (i < 1 || i > x@Dim[1] || j < 1 || j > x@Dim[2]) {
#               stop(gettextf("(i,j) should be (1..Dim[1],1..Dim[2]): i=%d, j=%d, Dim=(%d,%d)", i, j, x@Dim[1], x@Dim[2]))
#             }
#             u <- which(x@i == i & x@j ==j)
#             if (length(u) == 0) {
#               zeroM(x@si[x@i[u]], x@sj[x@j[u]])
#             } else if (length(u) == 1) {
#               x@x[[u]]
#             } else {
#               stop("Error")
#             }
#           })
#
# ## methods for [<-
#
# setReplaceMethod("[", signature(x = "bSparseMatrix", i = "index", j = "index", value = "Matrix"),
#                  function (x, i, j, ..., value) {
#                    stopifnot(length(i) == 1 && length(j) == 1)
#                    if (i < 1 || i > x@Dim[1] || j < 1 || j > x@Dim[2]) {
#                      stop(gettextf("(i,j) should be (1..Dim[1],1..Dim[2]): i=%d, j=%d, Dim=(%d,%d)", i, j, x@Dim[1], x@Dim[2]))
#                    }
#                    if (dim(value)[1] != x@si[i] || dim(value)[2] != x@sj[j]) {
#                      stop(gettextf("The (%d,%d)-element should be a %d-by-%d matrix. But the sieze of x is (%d,%d).", i, j, x@si[i], x@sj[j], dim(value)[1], dim(value)[2]))
#                    }
#                    u <- which(x@i == i & x@j ==j)
#                    if (length(u) == 0) {
#                      x@i <- c(x@i, as.integer(i))
#                      x@j <- c(x@j, as.integer(j))
#                      x@x <- c(x@x, list(value))
#                    } else if (length(u) == 1) {
#                      x@x[[u]] <- value
#                    } else {
#                      stop("Error")
#                    }
#                    x
#                  })
#
# setReplaceMethod("[", signature(x = "bSparseMatrix", i = "index", j = "index", value = "matrix"),
#                  function (x, i, j, ..., value) {
#                    stopifnot(length(i) == 1 && length(j) == 1)
#                    if (i < 1 || i > x@Dim[1] || j < 1 || j > x@Dim[2]) {
#                      stop(gettextf("(i,j) should be (1..Dim[1],1..Dim[2]): i=%d, j=%d, Dim=(%d,%d)", i, j, x@Dim[1], x@Dim[2]))
#                    }
#                    if (dim(value)[1] != x@si[i] || dim(value)[2] != x@sj[j]) {
#                      stop(gettextf("The (%d,%d)-element should be a %d-by-%d matrix. But the sieze of x is (%d,%d).", i, j, x@si[i], x@sj[j], dim(value)[1], dim(value)[2]))
#                    }
#                    u <- which(x@i == i & x@j ==j)
#                    if (length(u) == 0) {
#                      x@i <- c(x@i, as.integer(i))
#                      x@j <- c(x@j, as.integer(j))
#                      x@x <- c(x@x, list(value))
#                    } else if (length(u) == 1) {
#                      x@x[[u]] <- value
#                    } else {
#                      stop("Error")
#                    }
#                    x
#                  })
#
# setReplaceMethod("[", signature(x = "bSparseMatrix", i = "index", j = "index", value = "sSparseMatrix"),
#                  function (x, i, j, ..., value) {
#                    stopifnot(length(i) == 1 && length(j) == 1)
#                    if (i < 1 || i > x@Dim[1] || j < 1 || j > x@Dim[2]) {
#                      stop(gettextf("(i,j) should be (1..Dim[1],1..Dim[2]): i=%d, j=%d, Dim=(%d,%d)", i, j, x@Dim[1], x@Dim[2]))
#                    }
#                    if (dim(value)[1] != x@si[i] || dim(value)[2] != x@sj[j]) {
#                      stop(gettextf("The (%d,%d)-element should be a %d-by-%d matrix. But the sieze of x is (%d,%d).", i, j, x@si[i], x@sj[j], dim(value)[1], dim(value)[2]))
#                    }
#                    u <- which(x@i == i & x@j ==j)
#                    if (length(u) == 0) {
#                      x@i <- c(x@i, as.integer(i))
#                      x@j <- c(x@j, as.integer(j))
#                      x@x <- c(x@x, list(value))
#                    } else if (length(u) == 1) {
#                      x@x[[u]] <- value
#                    } else {
#                      stop("Error")
#                    }
#                    x
#                  })
#
# setReplaceMethod("[", signature(x = "bSparseMatrix", i = "index", j = "index", value = "NULL"),
#                  function (x, i, j, ..., value) {
#                    stopifnot(length(i) == 1 && length(j) == 1)
#                    if (i < 1 || i > x@Dim[1] || j < 1 || j > x@Dim[2]) {
#                      stop(gettextf("(i,j) should be (1..Dim[1],1..Dim[2]): i=%d, j=%d, Dim=(%d,%d)", i, j, x@Dim[1], x@Dim[2]))
#                    }
#                    u <- which(x@i == i & x@j ==j)
#                    if (length(u) == 0) {
#                    } else if (length(u) == 1) {
#                      x@i <- x@i[-u]
#                      x@j <- x@j[-u]
#                      x@x[[u]] <- NULL
#                    } else {
#                      stop("Error")
#                    }
#                    x
#                  })
#
# setAs("bSparseMatrix", "TsparseMatrix", function(from, to) {
#   as(convBtoM(from), "TsparseMatrix")
# })
#
# setAs("bSparseMatrix", "CsparseMatrix", function(from, to) {
#   as(convBtoM(from), "CsparseMatrix")
# })
#
# setAs("bSparseMatrix", "RsparseMatrix", function(from, to) {
#   as(convBtoM(from), "RsparseMatrix")
# })
#
# setMethod("print", signature(x = "bSparseMatrix"), function(x, ...) {
#   mat <- matrix(".", x@Dim[1], x@Dim[2])
#   for (u in 1:length(x@x)) {
#     mat[x@i[u], x@j[u]] <- gettextf(fmt="(%d,%d)", dim(x@x[[u]])[1], dim(x@x[[u]])[2])
#   }
#   mat <- data.frame(mat)
#   rownames(mat) <- sapply(1:x@Dim[1], function(i) gettextf("%d(%d)", i, x@si[i]))
#   colnames(mat) <- sapply(1:x@Dim[2], function(j) gettextf("%d(%d)", j, x@sj[j]))
#   print(mat)
# })
#
# ###
#
# setClass("kSparseMatrix", representation(x="list", OP="character"), contains="sSparseMatrix")
#
# kronMatrix <- function(..., op = "*") {
#   x <- list(...)
#   mdims <- c(1,1)
#   for (m in x) {
#     stopifnot(op == "*" || (op == "+" && dim(m)[1] == dim(m)[2]))
#     mdims <- mdims * dim(m)
#   }
#   km <- new("kSparseMatrix")
#   km@x <- x
#   km@MDim <- as.integer(mdims)
#   km@OP <- op
#   km
# }
#
# convKtoM <- function(km) {
#   mat <- c(1)
#   for (m in km@x) {
#     mat <- kronecker(X=mat, Y=as(m, "TsparseMatrix"))
#   }
#   mat
# }
#
# convKPtoM <- function(km) {
#   stopifnot(km@MDim[1] == km@MDim[2])
#   mat <- Matrix(0, km@MDim[1], km@MDim[2])
#   left <- 1
#   right <- km@MDim[1]
#   for (m in km@x) {
#     right <- right / dim(m)[1]
#     leftI <- sparseMatrix(i=1:left, j=1:left, x=rep.int(x=1, times=left), giveCsparse=TRUE)
#     rightI <- sparseMatrix(i=1:right, j=1:right, x=rep.int(x=1, times=right))
#     mat <- mat + kronecker(X=kronecker(X=leftI, Y=as(m, "TsparseMatrix")), Y=rightI, giveCsparse=TRUE)
#     left <- left * dim(m)[1]
#   }
#   mat
# }
#
# setAs("kSparseMatrix", "TsparseMatrix", function(from, to) {
#   if (from@OP == "*") {
#     as(convKtoM(from), "TsparseMatrix")
#   } else if (from@OP == "+") {
#     as(convKPtoM(from), "TsparseMatrix")
#   } else {
#     stop("Slot OP is neither '*' nor '+'.")
#   }
# })
#
# setAs("kSparseMatrix", "CsparseMatrix", function(from, to) {
#   if (from@OP == "*") {
#     as(convKtoM(from), "CsparseMatrix")
#   } else if (from@OP == "+") {
#     as(convKPtoM(from), "CsparseMatrix")
#   } else {
#     stop("Slot OP is neither '*' nor '+'.")
#   }
# })
#
# setAs("kSparseMatrix", "RsparseMatrix", function(from, to) {
#   if (from@OP == "*") {
#     as(convKtoM(from), "RsparseMatrix")
#   } else if (from@OP == "+") {
#     as(convKPtoM(from), "RsparseMatrix")
#   } else {
#     stop("Slot OP is neither '*' nor '+'.")
#   }
# })
#
# ###
#
# setClass("fSparseMatrix", representation(x="list", OP="character"), contains="sSparseMatrix")
#
# sumMatrix <- function(..., op = "+") {
#   x <- list(...)
#   mdims <- dim(x[[1]])
#   for (m in x) {
#     stopifnot(mdims[1] == dim(m)[1] && mdims[2] == dim(m)[2])
#   }
#   fm <- new("fSparseMatrix")
#   fm@x <- x
#   fm@MDim <- as.integer(mdims)
#   fm@OP <- op
#   fm
# }
#
# convFPtoM <- function(fm) {
#   mat <- Matrix(0, fm@MDim[1], fm@MDim[2])
#   for (m in fm@x) {
#     mat <- mat + as(m, "TsparseMatrix")
#   }
#   mat
# }
#
# convFMtoM <- function(fm) {
#   mat <- Matrix(0, fm@MDim[1], fm@MDim[2])
#   for (m in fm@x) {
#     mat <- mat - as(m, "TsparseMatrix")
#   }
#   mat
# }
#
# setAs("fSparseMatrix", "TsparseMatrix", function(from, to) {
#   switch(from@OP,
#          "+" = as(convFPtoM(from), "TsparseMatrix"),
#          "-" = as(convFMtoM(from), "TsparseMatrix"),
#          stop("Slot OP is neither '+' nor '-'."))
# })
#
# setAs("fSparseMatrix", "CsparseMatrix", function(from, to) {
#   switch(from@OP,
#          "+" = as(convFPtoM(from), "CsparseMatrix"),
#          "-" = as(convFMtoM(from), "CsparseMatrix"),
#          stop("Slot OP is neither '+' nor '-'."))
# })
#
# setAs("fSparseMatrix", "RsparseMatrix", function(from, to) {
#   switch(from@OP,
#          "+" = as(convFPtoM(from), "RsparseMatrix"),
#          "-" = as(convFMtoM(from), "RsparseMatrix"),
#          stop("Slot OP is neither '+' nor '-'."))
# })
