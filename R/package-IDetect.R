#' IDetect: Multiple generalised change-point detection using the Isolate-Detect methodology
#'
#' The \code{IDetect} package implements the Isolate-Detect methodology for
#' multiple generalised change-point detection, or sequence segmentation, in one-dimensional data
#' following the ``deterministic signal + noise'' model. The different structures that
#' are implemented are: piecewise-constant signal with Gaussian noise, piecewise-constant signal with
#' heavy tailed noise, piecewise-linear and continuous signal with Gaussian noise,
#' and piecewise-linear and continuous signal with heavy-tailed noise. The main routine
#' of the package is \code{\link{ID}}.
#'
#' @author Andreas Anastasiou, \email{a.anastasiou@lse.ac.uk}, Piotr Fryzlewicz, \email{p.fryzlewicz@lse.ac.uk}
#' @references ``Detecting multiple generalized change-points by isolating single ones'', Anastasiou
#'   and Fryzlewicz (2018), preprint.
#' @seealso
#' \code{\link{ID}}, \code{\link{ID_pcm}}, \code{\link{ID_cplm}}, \code{\link{ht_ID_pcm}},
#'   and \code{\link{ht_ID_cplm}}.
#' @examples
#' #See Examples for ID.
#' @docType package
#' @name IDetect
NULL
##NULL
