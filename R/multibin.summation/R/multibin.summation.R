#' Implements the proximal multibin summation (ProMBS) algorithm.
#'
#' \tabular{ll}{
#' Package: \tab multibin.summation\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0-1\cr
#' Date: \tab 2012-08-20\cr
#' URL: \tab http://www.adaptivesampling.org\cr
#' License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Inference in product partition models requires an average over
#' all possible partitions of the data. This algorithm focuses on
#' a particular case of product partition models where the data allows
#' a well-ordering for which it efficiently evaluates all possible
#' partitions.
#' 
#' The background and the algorithm are described in
#' Poppe, S, Benner, P, Elze, T. 
#'   A predictive approach to nonparametric inference for adaptive 
#'   sequential sampling of psychophysical experiments.
#'   Journal of Mathematical Psychology 56 (2012) 179-195
#'
#' See also
#' Daniel Barry and J. A. Hartigan
#'   Product Partition Models for Change Point Problems
#'   Ann. Statist. Volume 20, Number 1 (1992), 260-279.
#'
#' @name multibin.summation
#' @docType package
#' @title Implements the proximal multibin summation (ProMBS) algorithm.
#' @author Philipp Benner \email{Philipp.Benner@@mis.mpg.de}, Tobias Elze \email{Tobias.Elze@@schepens.harvard.edu}
#' @references
#'  Poppe, S, Benner, P, Elze, T. 
#'  A predictive approach to nonparametric inference for adaptive 
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @seealso \code{\link{prombs}} \code{\link{prombsExtended}}
#' @useDynLib multibin.summation
NULL
