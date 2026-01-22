#' Import STL triangles (facets) as Tibble
#'
#' Imports triangle centers and triangle normals of STL file as tibble. We use 
#' the word 'triangle' here to refer to the facets of an STL mesh to avoid 
#' confusion with the facets of compound eyes.
#'
#' @param file_name File name of STL to import.
#'
#' @return A tibble containing triangle centers and triangle normals of STL file.
#'
#' @export
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom dplyr filter pull select mutate arrange slice group_by ungroup left_join summarize distinct
#' first lead lag case_when bind_cols tibble as_tibble desc progress_estimated bind_rows all_of rename n 
#' mutate_all
#' @importFrom foreach foreach '%dopar%'
#' @importFrom magrittr '%>%'
#' @importFrom geometry dot
#' @importFrom graphics locator par abline hist
#' @importFrom grDevices grey.colors
#' @importFrom readr write_csv
#' @importFrom reshape2 melt
#' @importFrom rgl plot3d spheres3d selectpoints3d points3d text3d
#' @importFrom tidyr separate
#' @importFrom stats dist hclust setNames cutree median
#' @importFrom dplyr row_number
#' @examples
#' xxx: add example
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
  lines_tbl <- as_tibble(lines)
  if(verbose == TRUE){
    cat("Converting to tibble with", nrow(lines_tbl), "lines ...\n")
  }
  
  # get coordinates of triangle vertices from lines_tbl
  if(verbose == TRUE){
    cat("Extracting vertex coordinates of triangles...\n")
  }
  
  vertex_coords_triangles <- lines_tbl %>% 
    # slice(1:10) %>% 
    filter(grepl(".*vertex", value)) %>%
    mutate(value = gsub("^.*vertex", "vertex", value)) %>% 
    separate(value, into = c("value", "x", "y","z"), 
             sep = " ") %>% 
    select(-value) %>% 
    mutate_all(as.numeric)
  
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
    mutate("ID" = IDs) %>% 
    group_by(ID) %>% 
    mutate(x = mean(x), 
           y = mean(y),
           z = mean(z)) %>% 
    distinct(ID, .keep_all = T) %>%
    select(ID, x, y, z) %>% 
    arrange(ID) 
  
  # get normals of triangles
  if(verbose == TRUE){
    cat("Extracting triangle normals...\n")
  }
  
  normals <-  lines_tbl %>% 
    filter(grepl("facet normal", value)) %>% 
    separate(value, into = c("value.1", "value.2", "norm.x", "norm.y", "norm.z"), sep = " ") %>% 
    select(-c(value.1, value.2)) %>% 
    mutate_all(as.numeric) %>% 
    mutate(ID = row_number()) %>% 
    select(ID, norm.x, norm.y, norm.z)
  
  # join triangle center coordinates and their normals
  tri_centers_normals <- left_join(vertex_coords_triangle_centers, normals, by = "ID") %>% 
    ungroup()
  
  if(verbose == TRUE){
    cat("Found", nrow(tri_centers_normals), "triangle coordinates.\n")
  }
  
  if(plot_results == TRUE){
    # # plot triangle center coordinates
    plot3d(tri_centers_normals %>%
             select(x,y,z), 
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
