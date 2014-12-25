#' Load an R object file and return a specific object (or a list of all objects if name == NA)
#' without affecting the current environment.
#'
#' @param f file name
#' @param name object name to retrieve
load.objects <- function(f, name=NA) {
    e <- baseenv()
    if (is.na(name)) {
        name <- load(f, envir=e)
    }
    else {
        load(f, envir=e)
    }
    return(get(name, envir=e))
}

.wrap <- function(x) {
    s <- strsplit(x, "[ \t\n]", perl = TRUE, useBytes = TRUE)[[1]]
    w <- sapply(s, nchar) > 0
    paste(s[w], collapse=' ')
}

.my.dcast <- function(...) {
    d <- reshape2::dcast(...)
    rownames(d) <- d[,1]
    d[,2:ncol(d)]
}

# Merge two argument lists, replacing values in the first with
# those of the same name in the second, and adding to the first
# args that only appear in the second.
.merge.args <- function (x, y) {
    stopifnot(is.list(x), is.list(y))
    n <- names(y)
    n <- n[n != ""]
    for (ni in n) {
        x[[ni]] <- y[[ni]]
    }
    x
}

# Given a vector or matrix m, optionally reduce (if m is a matrix and indicies is not NULL)
# and then if the result is a vector, force it to be a matrix with 1 row (if margin == 1) or
# 1 column (if margin == 2).
.force.matrix <- function(m, margin=1, indices=NULL, default.name="x") {
    d <- NULL
    if (!is.null(indices)) {
        if (margin == 1) {
            n <- rownames(m)[indices]
            n[is.null(n)] <- default.name
            d <- list(n, colnames(m))
            m <- m[indices,]
        }
        else {
            n <- colnames(m)[indices]
            n[is.null(n)] <- default.name
            d <- list(rownames(m), n)
            m <- m[,indices]
        }
    }
    if (is.vector(m)) {
        if (margin == 1) {
            if (is.null(d)) {
                d <- list(default.name, names(m))
            }
            m <- matrix(m, nrow=1, ncol=length(m), dimnames=d)
        }
        else {
            if (is.null(d)) {
                d <- list(names(m), default.name)
            }
            m <- matrix(m, nrow=length(m), ncol=1, dimnames=d)
        }
    }
    m
}

# Collapse overlapping ranges. Expects a matrix or data
# frame with at least two columns: one defining the start
# of the range and the other defning the end.
.compress.intervals <- function(x, start.col=2, end.col=3, val.cols=NULL) {
    na.rows <- c()
    for (i in 1:(nrow(x)-1)) {
        if (x[i, end.col] >= (x[i+1, start.col] - 1)) {
            if (is.null(val.cols) || x[i, val.cols] == x[i+1, val.cols]) {
                x[i+1, start.col] <- x[i, start.col]
                na.rows <- c(na.rows, i)
            }
        }
    }
    if (length(na.rows) > 0) {
        x <- x[-na.rows,]
    }
    return(x)
}

.is.num <- function(x) !(is.na(x) | is.nan(x) | is.infinite(x))
