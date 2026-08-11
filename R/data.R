#' Example Compound-Eye Surface Mesh
#'
#' A reduced contiguous subset of a compound-eye surface mesh from an adult
#' *Drosophila melanogaster*. The source micro-computed tomography data are
#' described by Schoborg et al. (2019). The dataset contains triangle centres
#' and triangle surface normals and is intended for examples involving local
#' surface-height calculations.
#'
#' @format A data frame with 600 rows and 7 variables:
#' \describe{
#'   \item{ID}{Sequential triangle identifier.}
#'   \item{x, y, z}{Three-dimensional triangle-centre coordinates in micrometres.}
#'   \item{norm.x, norm.y, norm.z}{Components of the triangle surface normal.}
#' }
#'
#' @source Schoborg TA, Smith SL, Smith LN, Morris HD, Rusan NM (2019).
#'   Micro-computed tomography as a platform for exploring Drosophila
#'   development. *Development* 146(23): dev176685.
#'   \doi{10.1242/dev.176685}
"cv3d_example_surface"


#' Example Thresholded Local-Height Points
#'
#' A reduced subset of thresholded local-height points from the same adult
#' *Drosophila melanogaster* compound eye represented by
#' `cv3d_example_surface`. The source micro-computed tomography data are
#' described by Schoborg et al. (2019). The dataset is intended for examples of
#' facet-candidate condensation.
#'
#' @format A data frame with 74 rows and 5 variables:
#' \describe{
#'   \item{source_index}{Sequential point identifier.}
#'   \item{x, y, z}{Three-dimensional point coordinates in micrometres.}
#'   \item{height_value}{Thresholded local surface-height value.}
#' }
#'
#' @source Schoborg TA, Smith SL, Smith LN, Morris HD, Rusan NM (2019).
#'   Micro-computed tomography as a platform for exploring Drosophila
#'   development. *Development* 146(23): dev176685.
#'   \doi{10.1242/dev.176685}
"cv3d_example_thresholded"


#' Example Checked Facet Centres
#'
#' A reduced contiguous subset of checked compound-eye facet centres from an
#' adult *Drosophila melanogaster*. The source micro-computed tomography data
#' are described by Schoborg et al. (2019). The dataset is intended for
#' examples of facet-neighbour, facet-size, facet-normal, and optical-property
#' calculations.
#'
#' @format A data frame with 61 rows and 4 variables:
#' \describe{
#'   \item{ID}{Sequential facet identifier.}
#'   \item{x, y, z}{Three-dimensional facet-centre coordinates in micrometres.}
#' }
#'
#' @source Schoborg TA, Smith SL, Smith LN, Morris HD, Rusan NM (2019).
#'   Micro-computed tomography as a platform for exploring Drosophila
#'   development. *Development* 146(23): dev176685.
#'   \doi{10.1242/dev.176685}
"cv3d_example_facets"
