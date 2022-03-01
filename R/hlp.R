Hpd <- function(Val, Pr, post.mu, Cover)
{
    plog = data.frame(Val, Pr)
    names(plog) <- c("Val", "Pr")
    idx.mu <- which.min(abs(plog$Val - post.mu))
    n <- length(plog$Val)
    i = 1
    repeat{
        idx.ub <- idx.mu 
        idx.lb <- idx.mu - i
        if( sum(plog$Pr[idx.lb:idx.ub]) >= Cover )  break;
        idx.ub <- idx.mu + i
        if( sum(plog$Pr[idx.lb:idx.ub]) >= Cover )  break;
        i = i + 1
        if(idx.ub >= n ||  idx.lb <= 1) break; 
    }
    if(idx.ub >  n) idx.ub <- n
    if(idx.lb <  1) idx.lb <- 1
    UB = plog$Val[idx.ub]
    LB = plog$Val[idx.lb]
    Sum = sum(subset(plog, Val >= LB & Val<= UB)$Pr)
    c(LB,UB,Sum)
}

#' prior to posterior
#'
#' get point-wise posterior densities, given  data points with point-wise prior
#' probabilities, and a likelihood function.
#'
#' @param p point-wise prior probabilities (def=1)
#' @param f density function (i.e., [dnorm()])
#' @param ... point values to plug in likelihood function [f()]
#' @return posterior point densities
#'
#' The function realize Baysian theory in a numerical manner.
#' 
#' The data points in [...]  and  point-wise probabilities in [p] represent the
#' prior. The values in [p] does not have to adds up to 1.
#'
#' If the  probabilties [p] are not  give, an equal probable  prior is assumed,
#' that is, p=1 for all points.
#'
#' The likelihood function f() can be a density like [dnorm()] which takes data
#' and parameters from [...].
#'
#' The returned point-wise  probabilities adds up to 1, together  with the same
#' data points in [...], represent a posterior density.
p2p <- function(f, p=1, ...)
{
    p <- p * f(...)
    p / sum(p)
}

#' confidence peaks
#'
#' @param p point probabilities
#' @param c target confidence in [0, 1] (def = 0.95)
#' @param ... point values
#' @return mask of points with cumulative condidence just reaching \code{c}.
#'
#' @details the point probability and point values in (...) togetner approximate
#'     a density.
cfp <- function(p, c=.95, ...)
{
    p <- p / sum(p)
    m <- rep(TRUE, length(p))
    j <- order(p)
    m[seq_along(p)[j][cumsum(p[j]) < 1 - c]] <- FALSE
    m
}

#' confidence region
#'
#' @param r rally points, def = posterior mean
#' @param p point probabilities that adds up to 1;
#' @param c target confidence
#' @param ... point values
#' @return mask of points reaching condidence \code{c}.
cfr <- function(p, c=.95, r=pMu, ...)
{
    v <- cbind(...)
    if(is.null(r))
        r <- pMu
    if(is.character(r) || is.function(r))
        r <- do.call(r, list(p, ...))
    if(is.null(dim(r)))
        dim(r) <- c(length(r) / ncol(v), ncol(v))
    ## shortest distances from points to any rally.
    d <- euc(r, v)
    d <- apply(d, 2, min)
    ## mask
    m <- rep(TRUE, length(p))
    j <- order(d, decreasing = TRUE)
    m[seq_along(p)[j][cumsum(p[j]) < 1 - c]] <- FALSE
    m
}

#' posterior mean
#'
#' @param p point probabilities, adds up to 1.0.
pMu <- function(p, ...)
{
    colSums(p * cbind(...))
}

#' Euclidian distance (squared)
euc <- function(x, y)
{
    outer(rowSums(x^2), rowSums(y^2), `+`) - 2 * tcrossprod(x, y)
}
