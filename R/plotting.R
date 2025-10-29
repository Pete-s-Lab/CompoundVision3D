#' Plot eye in 3D
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
#' @examples
#' xxx: add example and change above parameter description
#' 
plot_eye <- function(facets_coords,
                     facet_colours = "red",
                     facet_sizes = 10,
                     facet_type = "p",
                     LM_coords = NULL,
                     LM_colours = "blue",
                     LM_sizes = 10,
                     LM_type = "p",
                     text_coords = NULL,
                     text_labels = NULL){
  
  # # testing
  # facets_coords = facets %>% 
  #   select(x, y, z)
  # facet_colours = facets$point_col
  # # facet_colours = "red"
  # facet_sizes = 60
  # facet_type = "p"
  # # facet_type = "s"
  # LM_coords = LMs %>% 
  #   select(x, y, z)
  # LM_colours = "blue"
  # # LM_colours = "red"
  # LM_sizes = 80
  # LM_type = "p"
  # # LM_type = "s"
  # text_coords = LMs %>% 
  #   select(x, y, z)
  # text_labels = LMs$ID
  
  # # close rgl window
  # while (rgl.cur() > 0) { close3d }
  
  if(facet_type == "s"){
    facet_radius = facet_sizes # /2
    plot3d(facets_coords, 
           col=facet_colours, 
           radius=facet_radius, 
           type = facet_type,
           aspect = "iso")
  } else if(facet_type == "p"){
    plot3d(facets_coords, 
           col=facet_colours, 
           size=facet_sizes, 
           type = facet_type,
           aspect = "iso")
  }
  
  if(!is.null(LM_coords)){
    if(LM_type == "s"){
      LM_radius = LM_sizes # /2
      spheres3d(LMs %>% 
                  select(x, y, z),
                col = "blue", radius=LM_radius)
    } else if(LM_type == "p"){
      points3d(LMs %>% 
                 select(x, y, z),
               col = "blue", size=LM_sizes)
    }
  }
  
  if(!is.null(LM_coords)){
    text3d(text_coords,
           texts = text_labels)
  }
}


#' create a continuous color ramp
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
#' @examples
#' xxx: add example and change above parameter description
#' 
continuous_color_ramp <- function(values, colors) {
  values[is.na(values)] <- min(values, na.rm = TRUE)
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}



#' Plot a facet with it's neighbours
#'
#' xxx
#'
#' @param df Dataframe with facet numbers and their coordinates.
#' @param facet The number of the facet.
#' @param radius Radius to plot the spheres.
#' @param threeD `Boolean` value for also plotting in 3D or just in 2D. Default: `FALSE` (= just 2D).
#'
#' @return xxx
#'
#' @export
#' @importFrom rgl plot3d spheres3d
#' @importFrom dplyr filter
#' @importFrom foreach foreach '%dopar%'
#' @examples
#' xxx: add example and change above parameter description
#' 
plot_facet_and_neighbours <- function(df,
                                      facet,
                                      radius = 1,
                                      threeD = FALSE){
  if(threeD == TRUE){
    require(rgl)
  }
  require(dplyr)
  
  # # testing
  # df = eye_L
  # facet = 1
  # radius = curr_facet_estimate/4
  
  facet_ <- facet
  curr_facet_coordinates <- df %>% 
    filter(facet == facet_)
  
  curr_neighbors <- as.numeric(str_split(df$neighbours[df$facet == facet], pattern = "; ")[[1]])
  curr_neighbors_coordinates <- df %>% 
    filter(facet %in% curr_neighbors)
  
  # identify the columns with the highest ranges for plotting
  column_ranges <- c(diff(range(curr_neighbors_coordinates$x)), diff(range(curr_neighbors_coordinates$y)), diff(range(curr_neighbors_coordinates$z)))
  names(column_ranges) <- c("x","y","z")
  x_y = column_ranges[1] * column_ranges[2]
  x_z = column_ranges[1] * column_ranges[3]
  y_z = column_ranges[2] * column_ranges[3]
  
  if(x_y > x_z & x_y > y_z){
    column1 <- names(column_ranges)[1]
    column2 <- names(column_ranges)[2]
  } else if(x_z > x_y & x_z > y_z){
    column1 <- names(column_ranges)[1]
    column2 <- names(column_ranges)[3]
  } else if(y_z > x_y & y_z > x_z){
    column1 <- names(column_ranges)[2]
    column2 <- names(column_ranges)[3]
  }
  
  plot(curr_facet_coordinates %>% 
         select(!!as.symbol(column1),
                !!as.symbol(column2)),
       cex = radius/2,
       col = "blue",
       pch=16,
  )
  points(curr_neighbors_coordinates %>% 
           select(!!as.symbol(column1),
                  !!as.symbol(column2)),
         cex = radius/2,
         pch=16,
         col = "orange")
  text(curr_facet_coordinates %>% 
         select(!!as.symbol(column1),
                !!as.symbol(column2)),
       labels = curr_facet_coordinates$facet,
       pos=1)
  text(curr_neighbors_coordinates %>% 
         select(!!as.symbol(column1),
                !!as.symbol(column2)),
       labels = curr_neighbors_coordinates$facet,
       pos=1)
  
  if(threeD == TRUE){
    close3d()
    plot3d(curr_facet_coordinates %>% 
             select(x,y,z),
           type = "s",
           radius = radius,
           col = "blue",
           aspect = "iso")
    spheres3d(curr_neighbors_coordinates %>% 
                select(x,y,z),
              radius = radius,
              col = "orange")
    text3d(curr_facet_coordinates %>% 
             select(x,y,z),
           texts = curr_facet_coordinates$facet,
           pos=1,
           cex=radius/2)
    text3d(curr_neighbors_coordinates %>% 
             select(x,y,z),
           texts = curr_neighbors_coordinates$facet,
           pos=1,
           cex=radius/2)
  }
}


#' Get optimal global viewing axis
#'
#' xxx
#'
#' @param df Dataframe with facet coordinates.
#' 
#' @return xxx
#'
#' @export
#' @importFrom dplyr pull
#' @examples
#' xxx: add example and change above parameter description
#' 
get_plot_axes <- function(df){
  # # testing
  # df = eye_L
  
  require(dplyr)
  
  range_x_y <- abs(diff(range(df %>% pull(x)))) * abs(diff(range(df %>% pull(y))))
  range_x_z <- abs(diff(range(df %>% pull(x)))) * abs(diff(range(df %>% pull(z))))
  range_y_z <- abs(diff(range(df %>% pull(y)))) * abs(diff(range(df %>% pull(z))))
  
  if(range_x_y > range_x_z & range_x_y > range_y_z){
    column1 = 1
    column2 = 2
    column3 = 3
  } else if(range_x_z > range_x_y & range_x_z > range_y_z){
    column1 = 1
    column2 = 3
    column3 = 2
  } else if(range_y_z > range_x_z & range_y_z > range_x_z){
    column1 = 2
    column2 = 3
    column3 = 1
  }
  
  return(c(column1, column2, column3))
}


#' Rotate the rgl viewport to face a specific direction vector
#'
#' xxx
#'
#' @param direction_vector Vector with three numbers indicating the desired 
#' viewing direction.
#' 
#' @return xxx
#'
#' @export
#' @importFrom dplyr pull
#' @examples
#' xxx: add example and change above parameter description
#' 
rotate_rgl_to_vector <- function(direction_vector) {
  # Normalize the direction vector
  direction_vector <- direction_vector / sqrt(sum(direction_vector^2))
  
  # Default "look-at" direction for the camera in rgl
  default_direction <- c(0, 0, -1)
  
  # Compute the rotation axis (cross product)
  rotation_axis <- c(
    default_direction[2] * direction_vector[3] - default_direction[3] * direction_vector[2],
    default_direction[3] * direction_vector[1] - default_direction[1] * direction_vector[3],
    default_direction[1] * direction_vector[2] - default_direction[2] * direction_vector[1]
  )
  
  # Compute the angle between the default and target directions (dot product)
  dot_product <- sum(default_direction * direction_vector)
  angle <- acos(dot_product)  # Angle in radians
  
  # Handle edge cases (parallel or antiparallel vectors)
  if (sum(rotation_axis^2) < 1e-6) {  # Nearly zero cross product -> no rotation needed
    if (dot_product < 0) {
      # Antiparallel case: Rotate 180 degrees around any perpendicular axis
      rotation_axis <- c(1, 0, 0)  # Choose an arbitrary perpendicular axis
      angle <- pi
    } else {
      # Parallel case: No rotation needed
      return()
    }
  }
  
  # Normalize the rotation axis
  rotation_axis <- rotation_axis / sqrt(sum(rotation_axis^2))
  
  # Rodrigues' rotation formula to compute the rotation matrix
  ux <- rotation_axis[1]
  uy <- rotation_axis[2]
  uz <- rotation_axis[3]
  cos_a <- cos(angle)
  sin_a <- sin(angle)
  R <- matrix(c(
    cos_a + ux^2 * (1 - cos_a), uy * ux * (1 - cos_a) + uz * sin_a, uz * ux * (1 - cos_a) - uy * sin_a,
    ux * uy * (1 - cos_a) - uz * sin_a, cos_a + uy^2 * (1 - cos_a), uz * uy * (1 - cos_a) + ux * sin_a,
    ux * uz * (1 - cos_a) + uy * sin_a, uy * uz * (1 - cos_a) - ux * sin_a, cos_a + uz^2 * (1 - cos_a)
  ), nrow = 3, byrow = TRUE)
  
  # Build the 4x4 transformation matrix
  userMatrix <- diag(4)
  userMatrix[1:3, 1:3] <- R
  
  # Apply the rotation
  par3d(userMatrix = userMatrix)
}

#' Calculate and plot facet vectors
#'
#' xxx
#'
#' @param df A `tibble` with the data columns (`x,y,z`) and the normal columns 
#' (`norm.x,norm.y,norm.z`)
#' @param vector_length_multipler a `numeric` value to change the length of the 
#' plotted vecots. Default: `3`.
#' 
#' @return xxx
#'
#' @export
#' @examples
#' xxx: add example and change above parameter description
#' 
make_segments <- function(df,
                          vector_length_multipler = 3){
  x_start <- df %>% pull(x)
  y_start <- df %>% pull(y)
  z_start <- df %>% pull(z)
  
  x_end <- df %>% pull(x) + df %>%
    pull(norm.x)*vector_length_multipler*(mean(df %>% pull(size), na.rm = TRUE))
  y_end <- df %>% pull(y) + df %>%
    pull(norm.y)*vector_length_multipler*(mean(df %>% pull(size), na.rm = TRUE))
  z_end <- df %>% pull(z) + df %>%
    pull(norm.z)*vector_length_multipler*(mean(df %>% pull(size), na.rm = TRUE))
  
  # Combine start and end points into a single vector
  segments <- cbind(x_start, y_start, z_start, x_end, y_end, z_end)
  
  segments3d(x=as.vector(t(segments[,c(1,4)])),
             y=as.vector(t(segments[,c(2,5)])),
             z=as.vector(t(segments[,c(3,6)])),
             col = viridis(n=6)[4],
             lwd=5)
  
  segments
}