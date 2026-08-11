#' Import Triangle Centres and Normals from an ASCII STL Mesh
#'
#' Reads an ASCII STL mesh and extracts the centre coordinates and surface
#' normal of each mesh triangle. One row is returned for each triangle in the
#' STL file. CV3D assumes that STL coordinate values are expressed in
#' micrometres (µm); the STL format itself does not encode a spatial unit.
#'
#' In CV3D, the term "triangle" refers to a triangular element of the STL mesh,
#' to distinguish these elements from the biological facets of a compound eye.
#'
#' @param file_name Path to an ASCII STL file.
#' @param plot_results Logical. If `TRUE`, plot the extracted triangle centres
#'   in an interactive 3D `rgl` window. Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress information during import.
#'   Default: `FALSE`.
#'
#' @return A tibble with one row per STL triangle and the columns: `ID`, a
#'   sequential triangle identifier; `x`, `y`, and `z`, the coordinates of the
#'   triangle centre; and `norm.x`, `norm.y`, and `norm.z`, the components of
#'   the triangle normal.
#'
#' @export
#' @examples
#' stl_file <- tempfile(fileext = ".stl")
#' writeLines(
#'   c(
#'     "solid example",
#'     "facet normal 0 0 1",
#'     "outer loop",
#'     "vertex 0 0 0",
#'     "vertex 1 0 0",
#'     "vertex 0 1 0",
#'     "endloop",
#'     "endfacet",
#'     "endsolid example"
#'   ),
#'   stl_file
#' )
#' triangles <- STL_triangles(stl_file)
#' triangles
#' unlink(stl_file)
#'
STL_triangles <- function(file_name, 
                          plot_results = FALSE,
                          verbose = FALSE){
  # dplyr NULLs
  ID <- x <- y <- z <- value <- value.1 <- value.2 <- row_number <- norm.x <- 
    norm.y <- norm.z <- NULL
  
  # # testing:
  # file_name <- "X:/Pub/2021/_Ruehr_AntVision/data/2_STLs/1_new/AV00001_Camponotus_hyatti_eye1_surface.stl"
  
  # load STL file as lines
  if(verbose == TRUE){
    cat("Importing", file_name, "...\n")
  }
  file_in <- file(file_name, open = "r")
  lines <- readLines(file_in)
  # delete first and last lines
  lines <- lines[-c(1, length(lines))]
  close(file_in)
  
  # convert character vector lines to tibble
  lines_tbl <- tibble::as_tibble(lines)
  if(verbose == TRUE){
    cat("Converting to tibble with", nrow(lines_tbl), "lines ...\n")
  }
  
  # get coordinates of triangle vertices from lines_tbl
  if(verbose == TRUE){
    cat("Extracting vertex coordinates of triangles...\n")
  }
  
  vertex_coords_triangles <- lines_tbl %>% 
    # slice(1:10) %>% 
    dplyr::filter(grepl(".*vertex", value)) %>%
    dplyr::mutate(value = gsub("^.*vertex", "vertex", value)) %>% 
    tidyr::separate(value, into = c("value", "x", "y","z"), 
             sep = " ") %>% 
    dplyr::select(-value) %>% 
    dplyr::mutate_all(as.numeric)
  
  # create vector of IDs (3 x same number per vertex coordinate)
  IDs <- c()
  for(i in 1:(nrow(vertex_coords_triangles)/3)){
    IDs[(length(IDs)+1):(length(IDs)+3)] <- rep(i, 3)
  }
  
  # get vertex coordinates of triangle centers
  if(verbose == TRUE){
    cat("Extracting coordinates of triangle centers...\n")
  }
  vertex_coords_triangle_centers <- vertex_coords_triangles %>% 
    dplyr::mutate("ID" = IDs) %>% 
    dplyr::group_by(ID) %>% 
    dplyr::mutate(x = mean(x), 
           y = mean(y),
           z = mean(z)) %>% 
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, x, y, z) %>% 
    dplyr::arrange(ID) 
  
  # get normals of triangles
  if(verbose == TRUE){
    cat("Extracting triangle normals...\n")
  }
  
  normals <-  lines_tbl %>% 
    dplyr::filter(grepl("facet normal", value)) %>% 
    tidyr::separate(value, into = c("value.1", "value.2", "norm.x", "norm.y", "norm.z"), sep = " ") %>% 
    dplyr::select(-c(value.1, value.2)) %>% 
    dplyr::mutate_all(as.numeric) %>% 
    dplyr::mutate(ID = dplyr::row_number()) %>% 
    dplyr::select(ID, norm.x, norm.y, norm.z)
  
  # join triangle center coordinates and their normals
  tri_centers_normals <- dplyr::left_join(vertex_coords_triangle_centers, normals, by = "ID") %>% 
    dplyr::ungroup()
  
  if(verbose == TRUE){
    cat("Found", nrow(tri_centers_normals), "triangle coordinates.\n")
  }
  
  if(plot_results == TRUE){
    # # plot triangle center coordinates
    if (!requireNamespace("rgl", quietly = TRUE)) stop("Package 'rgl' is required when plot_results = TRUE.", call. = FALSE)
    rgl::plot3d(tri_centers_normals %>%
             dplyr::select(x,y,z), 
           aspect = "iso", 
           col = "blue")
  }
  
  # # draw vectors
  # vec.mult <- 0.1
  # for(curr_facet in round(seq(1, nrow(tri_centers_normals), length.out = 150))){ # nrow(curr_facets)
  #   normal_vectors_df_subset <- tri_centers_normals %>%
  #     filter(ID == curr_facet) %>%
  #     select(norm.x, norm.y, norm.z)
  #   curr_facet_coordinates <- tri_centers_normals %>%
  #     filter(ID==curr_facet) %>%
  #     select(x,y,z)
  # 
  #   # find mean point of normalized normal vector ends
  #   norm.x <- normal_vectors_df_subset$norm.x
  #   norm.y <- normal_vectors_df_subset$norm.y
  #   norm.z <- normal_vectors_df_subset$norm.z
  # 
  #   lines3d(x = c(curr_facet_coordinates %>% pull(x), curr_facet_coordinates %>% pull(x) + norm.x*vec.mult),
  #           y = c(curr_facet_coordinates %>% pull(y), curr_facet_coordinates %>% pull(y) + norm.y*vec.mult),
  #           z = c(curr_facet_coordinates %>% pull(z), curr_facet_coordinates %>% pull(z) + norm.z*vec.mult),
  #           col = "red")
  # }
  if(verbose == TRUE){
    cat("done!\n")
  }
  return(tri_centers_normals)
}
