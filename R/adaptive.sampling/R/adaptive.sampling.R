#' Implements a predictive approach to non-parametric inference for adaptive sampling.
#'
#' \tabular{ll}{
#' Package: \tab adaptive.sampling\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0-1\cr
#' Date: \tab 2012-08-20\cr
#' URL: \tab http://www.adaptivesampling.org\cr
#' License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Implementing a predictive approach based on a hierarchical Bayesian
#' model, adaptive.sampling calculates utilities that can be used to
#' determine the successive step in an adaptive sequential sampling
#' measurement, such as the stimulus-response relation in a psychophysical
#' experiment.
#' 
#' The background and the algorithm are described in
#' Poppe, S, Benner, P, Elze, T. 
#'   A predictive approach to nonparametric inference for adaptive 
#'   sequential sampling of psychophysical experiments.
#'   Journal of Mathematical Psychology 56 (2012) 179-195
#'
#' @name adaptive.sampling
#' @docType package
#' @title Implements a predictive approach to non-parametric inference for adaptive sampling
#' @author Philipp Benner \email{Philipp.Benner@@mis.mpg.de}, Tobias Elze \email{Tobias.Elze@@schepens.harvard.edu}
#' @seealso \code{\link{adaptive.sampling}} \code{\link{binning.posterior}}
#' @references
#'  Poppe, S, Benner, P, Elze, T.
#'  A predictive approach to nonparametric inference for adaptive
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @useDynLib adaptive.sampling
NULL
