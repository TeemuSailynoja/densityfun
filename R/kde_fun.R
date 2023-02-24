#' Kernel density estimate function
#'
#' This function is used to obtain a callable KDE function, especially useful,
#' when the
#'
#' @param x the data from which the estimate is to be computed. For the default
#' method a numeric vector: long vectors are not supported.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such
#' that this is the standard deviation of the smoothing kernel.
#'
#' 'bw' can also be a character string giving a rule to choose the bandwidth.
#' See \link[stats]{bw.nrd}.
#'
#' The specified (or computed) value of 'bw' is multiplied by 'adjust'.
#' @param adjust the bandwidth used is actually adjust*bw. This makes it easy to
#'  specify values like ‘half the default’ bandwidth.
#' @param kernel a character string giving the smoothing kernel to be used. This
#'  must partially match one of "gaussian", "rectangular", "triangular",
#'  "epanechnikov", "biweight", "cosine" or "optcosine", with default
#'  "gaussian", and may be abbreviated to a unique prefix (single letter).
#'
#' "cosine" is smoother than "optcosine", which is the usual ‘cosine’ kernel in
#' the literature and almost MSE-efficient.
#' @param weights numeric vector of non-negative observation weights, hence of
#' same length as 'x'. The default 'NULL' is equivalent to
#' 'weights = rep(1/nx, nx)' where 'nx' is the length of (the finite entries of)
#'  'x[]'. If 'na.rm = TRUE' and there are NA's in 'x', they and the
#'  corresponding weights are removed before computations. In that case, when
#'  the original weights have summed to one, they are re-scaled to keep doing so.
#' @param na.rm logical; if TRUE, missing values are removed from x. If FALSE
#' any missing values cause an error.
#' @param ... further arguments for (non-default) methods.
#' @param fun.type Should a pdf or cdf be returned.
#'
#' @return a function with formal argument y specifying one or more evaluation
#' points for the pdf of cdf of the computed KDE.
#' @export
#'
kdefun <- function(x,
                    bw = "SJ",
                    adjust = 1,
                    kernel = c(
                      "gaussian",
                      "epanechnikov",
                      "rectangular",
                      "triangular",
                      "biweight",
                      "cosine",
                      "optcosine"
                    ),
                    weights = NULL,
                    na.rm = FALSE,
                    fun.type = c("pdf", "cdf"),
                    ...)
{
  chkDots(...)
  kernel <- match.arg(kernel)
  fun.type <- match.arg(fun.type)
  if (fun.type == "cdf")
    match.arg(kernel, c("gaussian"))
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  x <- as.vector(x)
  N <- length(x)
  if (has.wts <- !is.null(weights)) {
    if (length(weights) != N)
      stop("'x' and 'weights' have unequal length")
    some.na <- FALSE
  }
  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm) {
      N <- length(x <- x[!x.na])
      if (has.wts) {
        some.na <- TRUE
        trueD <- isTRUE(all.equal(1, sum(weights)))
        weights <- weights[!x.na]
        if (trueD)
          weights <- weights / sum(weights)
      }
    }
    else
      stop("'x' contains missing values")
  }
  nx <- N <- as.integer(N)
  if (is.na(N))
    stop(gettextf("invalid value of %s", "length(x)"), domain = NA)
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
    nx <- length(x)
  }
  if (!has.wts) {
    weights <- rep.int(1 / nx, nx)
    totMass <- nx / N
  }
  else {
    if (!all(is.finite(weights)))
      stop("'weights' must all be finite")
    if (any(weights < 0))
      stop("'weights' must not be negative")
    wsum <- sum(weights)
    if (any(!x.finite)) {
      weights <- weights[x.finite]
      totMass <- sum(weights) / wsum
    }
    else
      totMass <- 1
    if (!isTRUE(all.equal(1, wsum)))
      warning("sum(weights) != 1  -- will not get true density")
  }
  if (is.character(bw)) {
    if (nx < 2)
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(
      tolower(bw),
      nrd0 = stats::bw.nrd0(x),
      nrd = stats::bw.nrd(x),
      ucv = stats::bw.ucv(x),
      bcv = stats::bw.bcv(x),
      sj = ,
      `sj-ste` = stats::bw.SJ(x, method = "ste"),
      `sj-dpi` = stats::bw.SJ(x, method = "dpi"),
      stop("unknown bandwidth rule")
    )
  }
  if (!all(is.finite(bw)))
    stop("non-finite 'bw'")
  bw <- adjust * bw
  if (all(bw <= 0))
    stop("'bw' is not positive.")

  function(y) {
    switch(fun.type,
           pdf = sapply(y, function(y_i) {
             (1 / nx) * sum(switch(
               kernel,
               gaussian = stats::dnorm(y_i - x, sd = bw),
               rectangular = {
                 a <- bw * sqrt(3)
                 ifelse(abs(y_i - x) < a, 0.5 / a, 0)
               },
               triangular = {
                 a <- bw * sqrt(6)
                 ax <- abs(y_i - x)
                 ifelse(ax < a, (1 - ax / a) / a, 0)
               },
               epanechnikov = {
                 a <- bw * sqrt(5)
                 ax <- abs(y_i - x)
                 ifelse(ax < a, 3 / 4 * (1 - (ax / a) ^ 2) / a, 0)
               },
               biweight = {
                 a <- bw * sqrt(7)
                 ax <- abs(y_i - x)
                 ifelse(ax < a, 15 / 16 * (1 - (ax / a) ^ 2) ^ 2 / a, 0)
               },
               cosine = {
                 a <- bw / sqrt(1 / 3 - 2 / pi ^ 2)
                 ifelse(abs(y_i - x) < a, (1 + cos(pi * y_i - x / a)) /
                          (2 * a), 0)
               },
               optcosine = {
                 a <- bw / sqrt(1 - 8 / pi ^ 2)
                 ifelse(abs(y_i - x) < a,
                        pi / 4 * cos(pi * y_i - x / (2 * a)) / a,
                        0)
               }
             ))
           }),
           cdf =  sapply(y, function(y_i) {
             (1 / nx) * sum(switch(kernel, gaussian = stats::pnorm(y_i - x, sd = bw)))
           }))
  }
}

#' @describeIn kdefun KDE with adaptive bandwidth. After computing \eqn{ f_0 }, a regular KDE
#'  with bandwith determined by 'bw', each value in `x`, \eqn{ x_i } obtains an individual
#'  bandwidth equal to \eqn{h_i = \sqrt{\frac{\bar{f_0})}{f_0(x_i)}}.}
#'
#' @export
kdefun_adaptive <- function(x,
                   bw = "SJ",
                   adjust = 1,
                   kernel = c(
                     "gaussian",
                     "epanechnikov",
                     "rectangular",
                     "triangular",
                     "biweight",
                     "cosine",
                     "optcosine"
                   ),
                   weights = NULL,
                   na.rm = FALSE,
                   fun.type = c("pdf", "cdf"),
                   ...) {
  kde_tilde <- kdefun(x, bw, adjust, kernel, weights, na.rm)(x)
  lambda <- sqrt(mean(kde_tilde) / kde_tilde)
  kdefun(x, bw, adjust = lambda, kernel, weights, na.rm, fun.type)
}
