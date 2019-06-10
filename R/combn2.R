#Generate All Combinations of n Elements, Taken m at a Time
combn2 <- function (x, m, FUN = NULL, simplify = TRUE, ...) 
{
  stopifnot(length(m) == 1L, is.numeric(m))
  if (m < 0) 
    stop("m < 0", domain = NA)
  if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == 
      x) 
    x <- seq_len(x)
  n <- length(x)
  if (n < m) 
    stop("n < m", domain = NA)
  x0 <- x
  if (simplify) {
    if (is.factor(x)) 
      x <- as.integer(x)
  }
  m <- as.integer(m)
  e <- 0
  h <- m
  a <- seq_len(m)
  nofun <- is.null(FUN)
  if (!nofun && !is.function(FUN)) 
    stop("'FUN' must be a function or NULL")
  len.r <- length(r <- if (nofun) x[a] else FUN(x[a], ...))
  count <- as.integer(round(choose(n, m)))
  if (simplify) {
    dim.use <- if (nofun) 
      c(m, count)
    else {
      d <- dim(r)
      if (length(d) > 1L) 
        c(d, count)
      else if (len.r > 1L) 
        c(len.r, count)
      else c(d, count)
    }
  }
  if (simplify) 
    out <- matrix(r, nrow = len.r, ncol = count)
  else {
    out <- vector("list", count)
    out[[1L]] <- r
  }
  if (m > 0) {
    i <- 2L
    nmmp1 <- n - m + 1L
    while (a[1L] != nmmp1) {
      if (e < n - h) {
        h <- 1L
        e <- a[m]
        j <- 1L
      }
      else {
        e <- a[m - h]
        h <- h + 1L
        j <- 1L:h
      }
      a[m - h + j] <- e + j
      r <- if (nofun) 
        x[a]
      else FUN(x[a], ...)
      if (simplify) 
        out[, i] <- r
      else out[[i]] <- r
      i <- i + 1L
    }
  }
  if (simplify) {
    if (is.factor(x0)) {
      levels(out) <- levels(x0)
      class(out) <- class(x0)
    }
    dim(out) <- dim.use
  }
  out
}

# system.time(combn2(letters[1:4], 2))
# system.time(combn(letters[1:4], 2))
