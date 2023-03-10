# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title tetrascatt_c
#'
#' @description Computes scattering from a volumetric mesh
#' efficientlty, it is an auxiliary function called by tetrascatt function.
#' @seealso \code{\link{tetrascatt}}
#' @param  cw sound speed in the water in m/s
#'
#' @param  g  density constrast value, i.e  g= rho1/rhow, where rho1 and rhow
#' are the density values of the stcatterer and the unbounded medium
#' respectively
#'
#' @param  h  density sound speed constrast value   i.e  h= c1/cw, where c1 is
#' the sound speed of the stcatterer
#'
#' @param  freq an array of frequencies, where the scattering is computed.
#'
#' @param  Ver  a matrix with the vertex of the tetrahedra, each vertex has to
#' have three coordinates.
#' @param  Tet a matrix containing the four index of each tetrahedron.
#' @param  kversor  three component vector that indicates the direction of the
#' incident plane wave.
#' @return  A complex number array which contains the backward
#' differential far-field scattering cross-section (f infinity)
#' values at each frequency.
#' @export
tetrascatt_c <- function(cw, g, h, freq, Tet, Ver, kversor) {
    .Call(`_tetrascatt_tetrascatt_c`, cw, g, h, freq, Tet, Ver, kversor)
}

