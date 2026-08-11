#' CV3D: Three-Dimensional Analysis of Compound-Eye Surface Geometry
#'
#' CV3D provides tools for extracting and analysing compound-eye surface
#' geometry from three-dimensional data. The package implements the R-based
#' analysis stages of the CV3D workflow, including local surface-height
#' calculation, facet-candidate detection, facet-neighbour and size estimation,
#' facet surface-normal reconstruction, calculation of local optical
#' properties, landmark-based coordinate alignment, and spherical projection
#' of viewing directions.
#'
#' The package is designed to be used together with the CV3D UI, but its
#' individual functions can also be used directly in R.
#'
#' CV3D currently assumes all mesh and point-cloud coordinates are expressed in
#' micrometres (µm). Length-derived outputs therefore use micrometres; Snyder's
#' eye parameter is expressed in µm·rad.
#'
#' @aliases CV3D CompoundVision3D
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @keywords internal
"_PACKAGE"
