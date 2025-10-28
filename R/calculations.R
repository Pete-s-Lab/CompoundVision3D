#' Calculates the angle between 2 vectors
#'
#' Rotates 3D point cloud according to one defined vector so that this vector
#' is aligned to one of the global coordinate system axes.
#'
#' @param df A tibble containing coordinates in columns `x, y, z`.
#' @param line_points A 2x3 tibble containing coordinates of line to align the 
#' point cloud to. Must contain one row per point and columns `x, y, z`.
#' @param axis A character string defining the global axis to align to. Must be 
#' `x`, `y`, or `z`.
#' @return Returns a tibble with the aligned coordinates in columns `x, y, z`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
# calculate angle between 2 vectors (types: d = degrees, r = radians, both = c(regrees, radians))
calc_delta.phi <- function(curr.lens, curr.partner, type){
  #message('- ', curr.lens)
  #message('p = ', curr.partner)
  curr.x_norm_lens <- curr.lens[1]
  curr.y_norm_lens <- curr.lens[2]
  curr.z_norm_lens <- curr.lens[3]
  
  curr.x_norm_partner <- curr.partner[1]
  curr.y_norm_partner <- curr.partner[2]
  curr.z_norm_partner <- curr.partner[3]
  
  curr_delta_phi.rad <- acos(((curr.x_norm_lens * curr.x_norm_partner)
                              + (curr.y_norm_lens * curr.y_norm_partner)
                              + (curr.z_norm_lens * curr.z_norm_partner))
                             /
                               (sqrt(curr.x_norm_lens**2
                                     + curr.y_norm_lens**2
                                     + curr.z_norm_lens**2) *
                                  sqrt(curr.x_norm_partner**2
                                       + curr.y_norm_partner**2
                                       + curr.z_norm_partner**2)))
  curr_delta_phi.deg <- curr_delta_phi.rad*180/pi
  # message(paste0('-> ', curr_delta_phi.deg, '?'))
  if(type == "d"){
    return <- curr_delta_phi.deg
  } else if(type == "r"){
    return <- curr_delta_phi.rad
  } else if(type == "b"){
    return <- c(curr_delta_phi.deg, curr_delta_phi.rad)
  }
}


#' Finds the point closest to multiple lines
#'
#' Rotates 3D point cloud according to one defined vector so that this vector
#' is aligned to one of the global coordinate system axes.
#'
#' @param df A tibble containing coordinates in columns `x, y, z`.
#' @param line_points A 2x3 tibble containing coordinates of line to align the 
#' point cloud to. Must contain one row per point and columns `x, y, z`.
#' @param axis A character string defining the global axis to align to. Must be 
#' `x`, `y`, or `z`.
#' @return Returns a tibble with the aligned coordinates in columns `x, y, z`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
closest_point_to_lines <- function(lines) {
  require(MASS)
  n <- length(lines)
  
  # Construct the matrices A and b
  A <- matrix(0, nrow = 3, ncol = 3)
  b <- rep(0, 3)
  
  for (line in lines) {
    p <- line$point
    d <- line$direction
    
    # Project p onto the plane orthogonal to d
    P <- diag(3) - outer(d, d)
    A <- A + t(P) %*% P
    b <- b + t(P) %*% P %*% p
  }
  
  # Solve for the closest point
  closest_point <- solve(A, b)
  return(closest_point)
}

#' Translates coordinates according to a ROI
#'
#' Rotates 3D point cloud according to one defined vector so that this vector
#' is aligned to one of the global coordinate system axes.
#'
#' @param df A tibble containing coordinates in columns `x, y, z`.
#' @param line_points A 2x3 tibble containing coordinates of line to align the 
#' point cloud to. Must contain one row per point and columns `x, y, z`.
#' @param axis A character string defining the global axis to align to. Must be 
#' `x`, `y`, or `z`.
#' @return Returns a tibble with the aligned coordinates in columns `x, y, z`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
translate_ROIs <- function(df,
                           ROI_coordinates,
                           eye){
  # # testing
  # df = local_heights
  # # df = facet_positions_new
  # ROI_coordinates = curr_ROI_coordinates
  # eye = 2
  
  df_translated <- df %>% 
    mutate(x = x + (ROI_coordinates %>% 
                      filter(ROI == paste0("eye", eye)) %>% 
                      pull(x_coord) * px_size_eyes) -
             (ROI_coordinates %>% 
                filter(ROI == "head") %>% 
                pull(x_coord) * px_size_eyes),
           y = y + (ROI_coordinates %>% 
                      filter(ROI == paste0("eye", eye)) %>% 
                      pull(y_coord) * px_size_eyes) -
             (ROI_coordinates %>% 
                filter(ROI == "head") %>% 
                pull(y_coord) * px_size_eyes),
           z = z + (ROI_coordinates %>% 
                      filter(ROI == paste0("eye", eye)) %>% 
                      pull(z1_coord) * px_size_eyes) -
             (ROI_coordinates %>% 
                filter(ROI == "head") %>% 
                pull(z1_coord) * px_size_eyes))
  
  
  return(df_translated)
}

#' Extract ROI coordinates from crop_log file
#'
#' Rotates 3D point cloud according to one defined vector so that this vector
#' is aligned to one of the global coordinate system axes.
#'
#' @param df A tibble containing coordinates in columns `x, y, z`.
#' @param line_points A 2x3 tibble containing coordinates of line to align the 
#' point cloud to. Must contain one row per point and columns `x, y, z`.
#' @param axis A character string defining the global axis to align to. Must be 
#' `x`, `y`, or `z`.
#' @return Returns a tibble with the aligned coordinates in columns `x, y, z`.
#' @importFrom dplyr add_row
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
get_ROI_coordinates <- function(ROIs,
                                crop_log_data){
  results = tibble(ROI = character(),
                   x_coord = numeric(), y_coord = numeric(), 
                   x_length = numeric(), y_length = numeric(), 
                   z1_coord = numeric(), z2_coord = numeric())
  i=1
  for(i in 1:length(ROIs)){
    curr_ROI_name <- ROIs [i]
    curr_ROI <- crop_log_data$val[which(grepl(paste0("ROI_", curr_ROI_name), crop_log_data$var))]
    ROI_coords <- strsplit(gsub("\\);", "",
                                strsplit(curr_ROI, "\\(")[[1]][2]), ", ")[[1]]
    x_coord <- as.numeric(ROI_coords[1])
    y_coord <- as.numeric(ROI_coords[2])
    x_length <- as.numeric(ROI_coords[3])
    y_length <- as.numeric(ROI_coords[4])
    
    z1_coord <- as.numeric(crop_log_data$val[which(grepl(paste0("z_first_", curr_ROI_name), crop_log_data$var))])
    z2_coord <- as.numeric(crop_log_data$val[which(grepl(paste0("z_last_", curr_ROI_name), crop_log_data$var))])
    
    results = results %>%  
      add_row(ROI = curr_ROI_name,
              x_coord = x_coord, y_coord = y_coord, 
              x_length = x_length, y_length = y_length, 
              z1_coord = z1_coord, z2_coord = z2_coord)
  }
  return(results)
}




#' Calculate distance between two points in 3D.
#'
#' xxx: add description
#'
#' @param point1 A `vector` containing the `numeric` `x, y` and `z` coordinates 
#' pf point 1.
#' @param point2 A `vector` containing the `numeric` `x, y` and `z` coordinates 
#' pf point 2.
#' #' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`. 
#' @return Returns the `numeric` distance between point 1 and point 2.
#' @importFrom dplyr add_row
#'
#' @export
#' @examples
#' xxx: add example
#'
distance_3D <- function(point1, 
                        point2,
                        verbose = FALSE) {
  # Ensure the points are numeric vectors of length 3
  if (length(point1) != 3 || length(point2) != 3) {
    stop("Both points must be numeric vectors of length 3.")
  }
  
  # Calculate the differences in each dimension
  diff <- point2 - point1
  
  # Compute the squared differences and sum them up
  sum_of_squares <- sum(diff^2)
  
  # Take the square root of the sum of squared differences to get the distance
  distance <- sqrt(sum_of_squares)
  
  if(verbose == TRUE){
    cat("All done!\n")
  }
  return(distance)
}



#' Find facet neighbors and Calculate facet sizes
#'
#' xxx: add description
#'
#' @param df A `tibble` containing facet coordinates in columns `x, y, z`.
#' @param facet_size A `numeric` value containing the estimated facet size.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`. 
#' @return Returns a `tibble` containing the additional columns with info on 
#' facet size, facet neighbours and the number of neighbors of each facet.
#'
#' @export
#' @examples
#' xxx: add example
#'
find_neighbours <- function(df,
                            facet_size,
                            neighbour_threshold = 1.5,
                            cores = 1,
                            plot_file = NULL,
                            verbose = FALSE){
  require(doParallel)
  require(parallel)
  require(doSNOW)
  require(progress)
  
  # # testing
  # df = curr_facets
  # facet_size <- curr_facet_estimate
  # neighbour_threshold = neighbour_threshold
  # cores = 12
  # plot_file = file.path(facet_infos, 
  #                       gsub("_facets_LMs_combined.csv", "_neighbour_and_size_data.pdf", curr_file))
  # verbose = TRUE
  
  facet_size <- as.numeric(facet_size)
  cores <- as.numeric(cores)
  
  if(verbose == TRUE){
    cat("Creating distance matrix...\n")
  }
  facet_distance_matrix <- dist(df %>% 
                                  select(x,y,z),
                                method = "euclidean",
                                diag = TRUE)
  
  # transform distance matrix into data frame
  facet_distance_df_all <- suppressWarnings(melt(as.matrix((facet_distance_matrix), 
                                                           varnames = c("row", "col")))) %>% 
    as_tibble() %>% 
    filter(value > 0)
  
  colnames(facet_distance_df_all) <- c("facet_1", "facet_2", "distance")
  
  facet_distance_df <- facet_distance_df_all %>% 
    filter(distance <= facet_size * 3)
  
  # calculate facet sizes
  if(verbose == TRUE){
    cat(paste0("Calculating temporary sizes for ", nrow(df), " facets according to their closest 2 neighbours (multi-threaded)...\n"))
  }
  
  registerDoParallel(cores)
  facet = 50
  df_sizes <- foreach(facet = 1:nrow(df), # nrow(df)
                      .combine=rbind, .packages=c('dplyr')) %dopar% {
                        
                        facet_size <- facet_distance_df %>%
                          filter(facet_1 == df$facet[facet])  %>%
                          arrange(distance) %>%
                          slice(1:2) %>%
                          summarise(distance = mean(distance)) %>%
                          pull(distance)
                        
                        if(length(facet_size) == 0){
                          facet_size = 0
                        }
                        
                        # if (facet %% 10 == 0) {
                        #   cat(sprintf("Completed %d out of %d tasks\n", facet, nrow(df)))
                        # }
                        tmp <- facet_size
                        
                      }
  stopImplicitCluster()
  
  # remove facets where no size could have been calculated
  zero_sizes <- which(df_sizes==0)
  if(length(zero_sizes) > 0){
    warning("Distance data frame has ", length(which(df_sizes==0)), " zero size entries.")
  }
  
  # add size column to df
  df <- df %>% 
    mutate("size" = as.numeric(df_sizes))
  
  
  if(verbose == TRUE){
    cat("Finding six closest facets (multi-threaded)...\n")
  }
  
  facet=50
  registerDoParallel(cores)
  neighbour_columns <- foreach(facet = 1:nrow(df),
                               .combine=rbind, .packages=c('dplyr')) %dopar% {
                                 
                                 neighbour_list <- facet_distance_df %>% 
                                   filter(facet_1 == facet, facet_2 != facet) %>% 
                                   arrange(distance) %>% 
                                   slice_head(n=6) %>% 
                                   pull(facet_2) 
                                 
                                 tmp <- tibble(neighbours = paste(neighbour_list, collapse = "; "),
                                               number.of.neighbours = length(neighbour_list))
                               }
  stopImplicitCluster()
  
  
  df$neighbours <- neighbour_columns$neighbours
  # add a temprary 6, because all six closest facets have been taken so far
  df$number.of.neighbours <- neighbour_columns$number.of.neighbours
  
  
  # calculate median distance (~facet size) of all facets to their closest three neighbors xxx: this is median of all facets yet - bit maybe enough?
  median.size.tmp <- median(df$size, na.rm = TRUE)
  
  # find neighbors that are closer than median.size.tmp to current facet and remove them from facet.infos
  # but keep at least 3 neighbours
  if(verbose == TRUE){
    cat("Filtering close neighbours (multi-threaded)...\n")
  }
  
  
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  # progress bar ------------------------------------------------------------
  pb <- progress_bar$new(
    format = "facet = :facet [:bar] :elapsed | eta: :eta",
    total = nrow(df),    # 100 
    width = 60)
  
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(facet = df$facet[n]))
  } 
  
  opts <- list(progress = progress)
  l=100
  neighbours_tmp <- foreach(l = 1:nrow(df ),# nrow(df)
                            .combine=rbind,
                            .packages=c('dplyr'),
                            .options.snow = opts,
                            .errorhandling = "stop") %dopar% {
                              
                              
                              # registerDoParallel(cores)
                              # l=282
                              # neighbours_tmp <- foreach(l = 1:nrow(df), # nrow(df)
                              #                           .combine=rbind, .packages=c('dplyr')) %dopar% {
                              
                              curr_facet <- df$facet[l]
                              
                              neighbours_raw <- facet_distance_df %>%
                                filter(facet_1 == curr_facet, facet_2 != curr_facet) %>%
                                arrange(distance) %>%
                                filter(distance <= neighbour_threshold*median.size.tmp) %>% 
                                mutate(delta = distance - lag(distance)) 
                              
                              # increase neighbour_threshold until at least 2 neighbours were found
                              curr_neighbour_threshold <- neighbour_threshold
                              counter=0
                              while(nrow(neighbours_raw) < 2 & counter <= 50){
                                curr_neighbour_threshold <- 1.05*curr_neighbour_threshold
                                neighbours_raw <- facet_distance_df %>%
                                  filter(facet_1 == curr_facet, facet_2 != curr_facet) %>%
                                  arrange(distance) %>%
                                  filter(distance <= curr_neighbour_threshold*median.size.tmp) %>%
                                  mutate(delta = distance - lag(distance))
                                counter <- counter+1
                              }
                              
                              # # delta distances between neighbours
                              # neighbours_delta_threshold = neighbours_raw  %>% 
                              #   filter(delta != max(delta, na.rm = TRUE)) %>%
                              #   slice(-nrow(.)) %>% 
                              #   pull(delta) %>% 
                              #   median(na.rm=TRUE) * (2*neighbour_threshold)
                              
                              # maximum_allowed_row <- which(neighbours_raw$delta == max(neighbours_raw$delta, na.rm = TRUE))
                              # 
                              # neighbouring_facets <- neighbours_raw %>% 
                              #   slice(1:(maximum_allowed_row-1))
                              
                              # neighbouring_facets <- neighbours_raw %>%
                              #   filter(is.na(delta) | delta <= neighbours_delta_threshold) %>%
                              #   pull(facet_2)
                              
                              neighbouring_facets <- neighbours_raw %>%
                                pull(facet_2)
                              
                              # take a max. of 6 and min of 3 neighbors
                              if(length(neighbouring_facets) > 6){
                                neighbours_fin <- neighbouring_facets[1:6]
                                # } else if(length(neighbouring_facets) < 3){
                                # neighbours_fin <- facet_distance_df %>%
                                #   filter(facet_1 == curr_facet, facet_2 != curr_facet) %>%
                                #   arrange(distance) %>%
                                #   slice(1:3) %>%
                                #   pull(facet_2)
                              } else {
                                neighbours_fin <- neighbouring_facets
                              }
                              
                              tmp <- tibble(facet = curr_facet,
                                            neighbours = paste(neighbours_fin, collapse = "; "),
                                            number.of.neighbours = length(neighbours_fin))
                              # df$neighbours[df$facet == curr_facet] <- paste(neighbours, collapse = "; ")
                              # df$number.of.neighbours[df$facet == curr_facet] <- length(neighbours)
                            }
  # stopImplicitCluster()
  stopCluster(cl) 
  
  
  # add neighbours
  df <- df %>% 
    select(-c(neighbours, number.of.neighbours)) %>% 
    left_join(neighbours_tmp, by="facet")
  
  
  # # test plot
  # plot3d(df %>%
  #          select(x,y,z),
  #        radius = df$size/8,
  #        # size = mean(df$size)*1.2,
  #        col = viridis(n=6)[1],
  #        type="s",
  #        aspect = "iso")
  # 
  # texts3d(df %>%
  #           select(x,y,z),
  #         texts = df$facet,
  #         pos = 1)
  # tmp=1
  # for(tmp in 7:nrow(df)){
  #   plot_facet_and_neighbours(df,
  #                             facet = tmp,
  #                             radius = 2,
  #                             threeD = TRUE)
  # 
  #   invisible(readline(prompt="Press [enter] to continue"))
  # }
  
  
  # sort(df$number.of.neighbours)
  # hist(df$number.of.neighbours)
  
  # find facet sizes according to their neighbors
  if(verbose == TRUE){
    cat("Finding facet sizes according to closest neighbours (multi-threaded)...\n")
  }
  
  registerDoParallel(cores)
  u=1
  u=352
  eye_sizes <- foreach(u = 1:nrow(df), # nrow(df)
                       .combine=rbind, .packages=c('dplyr', 'filesstrings')) %dopar% {
                         curr_facet <- df$facet[u]
                         
                         curr.neighbors <- as.numeric(str_split(df$neighbours[df$facet == curr_facet], pattern = "; ")[[1]])
                         
                         # if(all(!is.na(curr.neighbors))){
                         curr_mean_distance <- facet_distance_df %>%
                           filter(facet_1 == curr_facet, facet_2 != curr_facet) %>%
                           filter(facet_2 %in% curr.neighbors) %>%
                           mutate(mean_distance = mean(distance)) %>%
                           pull(mean_distance) %>%
                           unique()
                         # } else{
                         #   curr_mean_distance <- 0
                         # }
                         
                         # df$size[df$facet == curr_facet] <- curr_mean_distance
                         tmp <- tibble(facet = curr_facet,
                                       size = curr_mean_distance)
                       }
  stopImplicitCluster()
  
  # hist(eye_sizes$size)
  
  # add sizes to df
  df <- df %>% 
    select(-size) %>% 
    left_join(eye_sizes %>% 
                select(facet, size),
              by="facet")
  
  # cleaning
  df <- df %>% 
    # select(-c(size)) %>% 
    # rename(size = size_final) %>%
    select(facet, neighbours, number.of.neighbours, size)
  
  # calculate mean of sizes according to facet neighbors
  if(verbose == TRUE){
    cat("Calculating mean sizes according to all their actual neighbors...\n")
  }
  
  df$mean_size <- NA
  for(l in 1:nrow(df)){
    curr_facet <- df$facet[l]
    curr_neighbours <- as.numeric(str_split(df$neighbours[df$facet == curr_facet], pattern = "; ")[[1]])
    
    # get their sizes
    curr_mean_size <- df %>% 
      filter(facet %in% curr_neighbours) %>% 
      summarise(mean_size = mean(size)) %>% 
      pull(mean_size)
    
    df$mean_size[l] <- curr_mean_size
    
  }
  
  # remove temporary size column
  df <- df %>% 
    select(-size) %>% 
    dplyr::rename(size = mean_size)
  
  
  if(verbose == TRUE){
    cat("Plotting infos to plot device...\n")
    # cat(min(df$size), "\n")
    # cat(max(df$size), "\n")
  }
  par(mfrow = c(2,1))
  hist(df$number.of.neighbours, 
       breaks = c(0:7),
       main = "Number of neighbours",
       xlab = "Number of neighbours")
  hist(df$size, 
       # breaks = seq(min(df$size), max(df$size), 
       #              length.out=16),
       main = "Raw facet sizes",
       xlab = "Facet size (um)")
  par(mfrow = c(1,1))
  
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("Plotting infos to", plot_file, "\n")
    }
    
    # PDF plots
    pdf(plot_file, # , today()
        onefile = TRUE, paper = "a4")
    
    par(mfrow = c(2,1))
    hist(df$number.of.neighbours, 
         breaks = c(0:7),
         main = "Number of neighbours",
         xlab = "Number of neighbours")
    hist(df$size, 
         # breaks = seq(min(df$size), max(df$size), 
         #              length.out=16),
         main = "Raw facet sizes",
         xlab = "Facet size (um)")
    par(mfrow = c(1,1))
    
    dev.off()
  }
  
  
  if(verbose == TRUE){
    cat("All done!\n")
  }
  
  return(df)
}



#' Calculate facet normals according to their spacial positions
#'
#' xxx: add description
#'
#' @param df A tibble containing facet coordinates in columns `x, y, z`.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`.
#' @return Returns a `tibble` containing the additional columns with info on 
#' facet normals in x, y, and z direction for each facet.
#'
#' @export
#' @examples
#' xxx: add example
#'
get_facet_normals <- function(df,
                              cores = 1,
                              plot_file = NULL,
                              verbose = FALSE){
  
  require(parallel)
  require(doSNOW)
  require(progress)
  
  # # testing
  # df = eye_L
  # cores = 12
  # verbose = TRUE
  # plot_file = file.path(facet_infos,
  #                       gsub("_facets_LMs_combined.csv", "_normal_data.pdf", curr_file))
  
  # get mean coordinate of facets
  eyes_mean <- df %>% 
    mutate(mean_x = mean(x),
           mean_y = mean(y),
           mean_z = mean(z)) %>%
    distinct(mean_x, mean_y, mean_z) 
  
  
  # plot3d(df %>%
  #          select(x,y,z),
  #        aspect = "iso")
  # text3d(df %>%
  #          select(x,y,z),
  #        texts = df %>%
  #          pull(facet))
  # points3d(eyes_mean, size=20)
  
  
  # remove neighbours that are not part of df anymore
  # iterate until no more rows need to be deleted
  # df_tmp <- df
  # df <- df_tmp
  curr_difference <- 1
  counter <- 0
  while(curr_difference != 0){
    counter <- counter+1
    # cat(counter, "\n)
    # for(m in 1:2){
    l=202
    facets_to_remove <- c()
    for(l in 1:nrow(df)){
      curr_facet <- df$facet[l]
      curr_neighbours <- as.numeric(str_split(df$neighbours[l], pattern = "; ")[[1]])
      
      neighbours_to_keep <- curr_neighbours[which(curr_neighbours %in% df$facet)]
      neighbours_to_remove <- curr_neighbours[which(curr_neighbours %in% df$facet == FALSE)]
      facets_to_remove <- c(facets_to_remove, neighbours_to_remove)
      # if(length(neighbours_to_remove > 0)) cat(neighbours_to_remove,"\n")
      
      df$neighbours[l] <- paste(neighbours_to_keep, collapse = "; ")
      df$number.of.neighbours[l] <- length(neighbours_to_keep)
    }
    
    last_nrow <- nrow(df)
    df <- df %>% 
      filter(number.of.neighbours>=2)
    
    curr_nrow <- nrow(df)
    curr_difference <- last_nrow-curr_nrow
  }
  
  # calculate facet normals according to their neighbours
  if(verbose == TRUE){
    cat(paste0("Calculating ", nrow(df), " facet normals according to their neighbors' coordinates (multi-threaded)...\n"))
    cat(paste0("This may take a while, because ", nrow(df), " x ", sum(df$number.of.neighbours), " = ", nrow(df)*sum(df$number.of.neighbours), " calculations will be performed.\n"))
  }
  
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  # progress bar ------------------------------------------------------------
  pb <- progress_bar$new(
    format = "facet = :facet [:bar] :elapsed | eta: :eta",
    total = nrow(df),    # 100 
    width = 60)
  
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick(tokens = list(facet = df$facet[n]))
  } 
  
  opts <- list(progress = progress)
  l=308
  df_normals <- foreach(l = 1:nrow(df),# nrow(df)
                        .combine=rbind, 
                        .packages=c('dplyr', 'stringr', 'CompoundVision3D'), 
                        .options.snow = opts,
                        .errorhandling = "stop") %dopar% {
                          
                          # print(l)
                          curr_facet <- df$facet[l]
                          curr_neighbours <- as.numeric(str_split(df$neighbours[l], pattern = "; ")[[1]])
                          
                          # if(all(!is.na(curr_neighbours))){
                          
                          # get coordinates of current facet and neighbors
                          coords_facet <- df %>% 
                            filter(facet == curr_facet) %>% 
                            select(x,y,z) %>% 
                            unlist()
                          
                          # calculate normal from eye center to current facet
                          center_vector <- coords_facet %>% unlist() - eyes_mean %>% unlist()
                          
                          # normalize center_vector
                          center_vector <- center_vector / sqrt(sum(center_vector^2))
                          # print(center_vector)
                          
                          
                          # # plot vector from center of eyes to curr_facet
                          # lines3d(x = c(eyes_mean %>% pull(mean_x), eyes_mean %>% pull(mean_x)+center_vector[1]*300),
                          #         y = c(eyes_mean %>% pull(mean_y), eyes_mean %>% pull(mean_y)+center_vector[2]*300),
                          #         z = c(eyes_mean %>% pull(mean_z), eyes_mean %>% pull(mean_z)+center_vector[3]*300),
                          #         col = "red")
                          
                          # get distances between curr_neighbours
                          curr_dists <- tibble(
                            n1 = numeric(),
                            n2 = numeric(),
                            dist = numeric()
                          )
                          n1=1
                          n2=2
                          for(n1 in 1:length(curr_neighbours)){
                            for(n2 in 1:length(curr_neighbours)){
                              curr_neighbour_1 <- curr_neighbours[n1]
                              curr_neighbour_2 <- curr_neighbours[n2]
                              coords_neighbour_1 <- df %>%
                                filter(facet == curr_neighbour_1) %>%
                                select(x,y,z) %>%
                                unlist()
                              coords_neighbour_2 <- df %>%
                                filter(facet == curr_neighbour_2) %>%
                                select(x,y,z) %>%
                                unlist()
                              curr_dist <- distance_3D(point1 = coords_neighbour_1,
                                                       point2 = coords_neighbour_2,
                                                       verbose = FALSE)
                              
                              curr_dists <- curr_dists %>% 
                                add_row(n1 = curr_neighbour_1,
                                        n2 = curr_neighbour_2,
                                        dist = curr_dist)
                            }
                          }
                          
                          # clean list from duplicates
                          curr_dists_clean <- curr_dists %>%
                            filter(dist > 0) %>% 
                            distinct(dist, .keep_all = TRUE) %>% 
                            arrange(dist) %>%
                            # mutate(delta_dist = dist - lag(dist, default = dist[1]))
                            slice(1:length(curr_neighbours))
                          
                          # print(curr_dists_clean)
                          
                          
                          # create triangle with curr_neighbor and get normal
                          curr_normals_x <- c()
                          curr_normals_y <- c()
                          curr_normals_z <- c()
                          curr_normals_angles <- c()
                          
                          n=1
                          for(n in 1:nrow(curr_dists_clean)){
                            # get neighbors for current triangle
                            curr_neighbour_1 <- curr_dists_clean %>%
                              slice(n) %>%
                              pull(n1)
                            curr_neighbour_2 <- curr_dists_clean %>%
                              slice(n) %>%
                              pull(n2)
                            
                            # print(paste0("Building triangle with facets ", curr_facet, ", ", curr_neighbour_1, " and ", curr_neighbour_2))
                            coords_neighbour_1 <- df %>%
                              filter(facet == curr_neighbour_1) %>%
                              select(x,y,z) %>%
                              unlist()
                            coords_neighbour_2 <- df %>%
                              filter(facet == curr_neighbour_2) %>%
                              select(x,y,z) %>%
                              unlist()
                            
                            curr_normal <- calculate_normal(A = coords_facet,
                                                            B = coords_neighbour_1,
                                                            C = coords_neighbour_2,
                                                            normalize = TRUE)
                            
                            # # check if facet normal points in same direction as triangle normal
                            # print(curr_normal)
                            # print(center_vector)
                            
                            curr_angle <- angle_between_vectors(a = curr_normal,
                                                                b = center_vector)
                            # print(curr_angle)
                            
                            # check if normals point to same direction
                            if(curr_angle < 0){
                              curr_normal <- -1*curr_normal
                            }
                            
                            # print(paste0("Adding normals from triangle with facets ", curr_facet, ", ", curr_neighbour_1, " and ", curr_neighbour_2))
                            curr_normals_x <- c(curr_normals_x, curr_normal[1])
                            curr_normals_y <- c(curr_normals_y, curr_normal[2])
                            curr_normals_z <- c(curr_normals_z, curr_normal[3])
                            curr_normals_angles <- c(curr_normals_angles, curr_angle)
                          } 
                          
                          tmp <- tibble(facet = curr_facet,
                                        norm.x = mean(curr_normals_x[1]), 
                                        norm.y = mean(curr_normals_y[1]), 
                                        norm.z = mean(curr_normals_z[1]))
                          # print("******************************")
                        }
  
  stopCluster(cl) 
  
  
  # average facet angles
  l=1
  df_normals_avg <- tibble(facet = numeric(),
                           norm.x = numeric(),
                           norm.y = numeric(),
                           norm.z = numeric())
  
  for(l in 1:nrow(df)) {
    
    # print(l)
    curr_facet <- df$facet[l]
    curr_neighbours <- as.numeric(str_split(df$neighbours[l], pattern = "; ")[[1]])
    
    curr_normal <- df_normals %>% 
      filter(facet == curr_facet) %>% 
      select(norm.x,norm.y,norm.z)
    
    curr_neighbor_normals <- df_normals %>% 
      filter(facet %in% curr_neighbours) %>% 
      select(norm.x,norm.y,norm.z)
    
    curr_normal_avg_x = mean(c(rep(curr_normal$norm.x, 1),
                               curr_neighbor_normals$norm.x))
    curr_normal_avg_y = mean(c(rep(curr_normal$norm.y, 1),
                               curr_neighbor_normals$norm.y))
    curr_normal_avg_z = mean(c(rep(curr_normal$norm.z, 1),
                               curr_neighbor_normals$norm.z))
    
    df_normals_avg <- df_normals_avg %>% 
      add_row(facet = curr_facet,
              norm.x = curr_normal_avg_x,
              norm.y = curr_normal_avg_y,
              norm.z = curr_normal_avg_z)
  }
  
  # replace original values
  df_normals <- df_normals_avg
  
  if(verbose == TRUE){
    cat("Plotting infos to plot device...\n")
  }
  par(mfrow = c(3,1))
  hist(df_normals$norm.x, 
       breaks = seq(min(df_normals$norm.x), max(df_normals$norm.x), 
                    length.out=16),
       main = "x",
       xlab = "normals x")
  hist(df_normals$norm.y, 
       breaks = seq(min(df_normals$norm.y), max(df_normals$norm.y), 
                    length.out=16),
       main = "y",
       xlab = "normals y")
  hist(df_normals$norm.z, 
       breaks = seq(min(df_normals$norm.z), max(df_normals$norm.z), 
                    length.out=16),
       main = "z",
       xlab = "normals z")
  par(mfrow = c(1,1))
  
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("Plotting infos to", plot_file, "\n")
    }
    
    # PDF plots
    pdf(plot_file, # , today()
        onefile = TRUE, paper = "a4")
    
    par(mfrow = c(3,1))
    hist(df_normals$norm.x, 
         breaks = seq(min(df_normals$norm.x), max(df_normals$norm.x), 
                      length.out=16),
         main = "x",
         xlab = "normals x")
    hist(df_normals$norm.y, 
         breaks = seq(min(df_normals$norm.y), max(df_normals$norm.y), 
                      length.out=16),
         main = "y",
         xlab = "normals y")
    hist(df_normals$norm.z, 
         breaks = seq(min(df_normals$norm.z), max(df_normals$norm.z), 
                      length.out=16),
         main = "z",
         xlab = "normals z")
    par(mfrow = c(1,1))
    
    dev.off()
  }
  
  if(verbose == TRUE){
    cat("All done!\n")
  }
  
  return(df_normals)
}



#' Get optic parameter approximations
#'
#' xxx: add description
#'
#' @param df A tibble containing facet coordinates in columns `x, y, z`.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`.
#' @return Returns a `tibble` containing the additional columns with info on 
#' the Eye Parameter (P), the inter-facet angle (delta.phi) and acuity (CPD) for
#' each facet.
#'
#' @export
#' @examples
#' xxx: add example
#'
get_optic_properties <- function(df,
                                 cores = 1,
                                 plot_file = NULL,
                                 verbose = FALSE){
  
  require(doParallel)
  
  # # testing
  # df = eye_L
  # cores = 12
  # plot_file = file.path(facet_infos,
  #                       gsub("_facets_LMs_combined.csv", "_optics_parameters.pdf", curr_file))
  # verbose = TRUE
  
  if(verbose == TRUE){
    cat(paste0("Calculating IF angles, P, v, and CPD for ", nrow(df), " facets (multi-threaded)...\n"))
  }
  
  # eye_center = colMeans(df %>% 
  #                         select(x,y,z), na.rm = TRUE)
  # c(df %>% pull(x) %>% mean(),
  #              df %>% pull(y) %>% mean(),
  #              df %>% pull(z) %>% mean())
  
  # mean_eye_normal =colMeans(df %>% 
  #                             select(norm.x,norm.y,norm.z), na.rm = TRUE)
  # c(df %>% pull(norm.x) %>% mean(., na.rm = TRUE),
  #                   df %>% pull(norm.y) %>% mean(., na.rm = TRUE),
  #                   df %>% pull(norm.z) %>% mean(., na.rm = TRUE))
  
  registerDoParallel(cores)
  l=1
  dphi_Ps_CPDs <- foreach(l = 1:nrow(df), # nrow(df)
                          .combine=rbind, .packages=c('dplyr', 'filesstrings', 'CompoundVision3D')) %dopar% {
                            
                            facet_no <- df$facet[l]
                            
                            curr_facet_coords <- df %>%
                              filter(facet == facet_no) %>%
                              select(x, y, z)
                            
                            curr_facet_normal <- df %>%
                              filter(facet == facet_no) %>%
                              select(norm.x, norm.y, norm.z)
                            if(all(!is.na(curr_facet_normal))){
                              if(all(curr_facet_normal != 0)){
                                # get neighbouring facets without NAs
                                curr_facet_neighbours <- as.numeric(
                                  str_split(df %>%
                                              filter(facet==facet_no) %>%
                                              pull(neighbours),
                                            pattern = "; ")[[1]])
                                
                                if(all(!is.na(curr_facet_neighbours))){
                                  
                                  delta_phis.rad <- c()
                                  delta_phis.deg <- c()
                                  
                                  neighbour <- curr_facet_neighbours[3]
                                  for (neighbour in curr_facet_neighbours) {
                                    curr_neighbour_coords <- df %>%
                                      filter(facet == neighbour) %>%
                                      select(x, y, z)
                                    
                                    curr_neighbour_normal <- df %>%
                                      filter(facet == neighbour) %>%
                                      select(norm.x, norm.y, norm.z)
                                    
                                    
                                    curr_delta_phi.rad <- calc_delta.phi(curr_facet_normal, curr_neighbour_normal, type = "r")
                                    if(is.na(curr_delta_phi.rad)){
                                      curr_delta_phi.rad = data.frame(0,0)
                                    }
                                    curr_delta_phi.deg <- curr_delta_phi.rad*180/pi
                                    
                                    
                                    delta_phis.rad <- c(delta_phis.rad, curr_delta_phi.rad[1,1])
                                    delta_phis.deg <- c(delta_phis.deg, curr_delta_phi.deg[1,1])
                                    
                                    # print(paste0("curr. ° = ", curr_delta_phi.deg))
                                  }
                                  
                                  # df$delta_phi.rad[l] <- mean(delta_phis.rad)
                                  # df$delta_phi.deg[l] <- mean(delta_phis.deg)
                                  
                                  curr_P <- df$size[as.numeric(as.character(facet_no))] * mean(delta_phis.rad)  * (sqrt(3)/2) # eye parameter Snyder 1977. Brigitte: (sqrt(3)/2) *
                                  # sampling frequency calculated after Feller et al. 2021 as CPD. Snyder 1977 for hexagonal lattice of visual axes: 1/sqrt(3) *  mean(delta_phis.rad)
                                  curr_CPD <- 1 / (2 * mean(delta_phis.rad))
                                } else {
                                  warning("No neighbor data for facet ", facet_no)
                                  curr_P <- 0
                                  curr_CPD <- 0
                                }
                              } else {
                                warning("No neighbor data for facet ", facet_no)
                                curr_P <- 0
                                curr_CPD <- 0
                              }
                            } else{
                              warning("No neighbor data for facet ", facet_no)
                              curr_P <- 0
                              curr_CPD <- 0
                            }
                            
                            tmp <- rbind(tibble(facet = facet_no, delta_phi.rad =  mean(delta_phis.rad), delta_phi.deg = mean(delta_phis.deg), P = curr_P, CPD = curr_CPD))
                            # df$P[l] <- curr_P
                            # df$CPD[l] <- curr_CPD
                          }
  stopImplicitCluster()
  
  # adding size and neighbour column to dphi_Ps_CPDs
  dphi_Ps_CPDs <- dphi_Ps_CPDs %>% 
    left_join(df %>% 
                select(facet, size, neighbours), 
              by = "facet")
  
  # adding neighbour column to dphi_Ps_CPDs
  dphi_Ps_CPDs <- dphi_Ps_CPDs %>% 
    left_join(df %>% 
                select(facet, neighbours))
  
  
  if(verbose == TRUE){
    cat("Calculating weighted mean IF angles...\n")
  }
  
  # find mean IF angle for each facet
  dphi_Ps_CPDs$mean.delta_phi.deg <- NA
  q = 1
  for(q in 1:nrow(dphi_Ps_CPDs)){
    curr_facet <- dphi_Ps_CPDs$facet[q]
    curr.neighbors <- as.numeric(
      str_split(dphi_Ps_CPDs %>% 
                  filter(facet==curr_facet) %>% 
                  pull(neighbours), 
                pattern = "; ")[[1]])
    
    # get number of neghbors for weighted averaging
    no_of_neihbors <- length(curr.neighbors)
    
    weighting_factor <- 6 - no_of_neihbors
    
    # get mean of IF-angles with curr facet angle doubled-weighted
    dphi_Ps_CPDs$mean.delta_phi.deg[q] <- mean(c(rep(dphi_Ps_CPDs$delta_phi.deg[curr_facet], (weighting_factor+1)), 
                                                 dphi_Ps_CPDs$delta_phi.deg[curr_facet],
                                                 dphi_Ps_CPDs$delta_phi.deg[curr.neighbors]))
  }
  dphi_Ps_CPDs$mean.delta_phi.rad <- dphi_Ps_CPDs$mean.delta_phi.deg*pi/180
  
  # calculate P and CPD with average angles per facet
  if(verbose == TRUE){
    cat("Calculating mean P and CPD angles...\n")
  }
  dphi_Ps_CPDs$P.mean <- NA
  dphi_Ps_CPDs$CPD.mean <- NA
  dphi_Ps_CPDs$v.mean <- NA
  facet = 1
  for(facet in 1:nrow(dphi_Ps_CPDs)){
    dphi_Ps_CPDs$P.mean[facet] <- dphi_Ps_CPDs$size[facet] * 
      dphi_Ps_CPDs$mean.delta_phi.rad[facet] * 
      (sqrt(3)/2) # eye parameter Snyder 1977. Brigitte: (sqrt(3)/2) * 
    dphi_Ps_CPDs$CPD.mean[facet] <- 1 / (2 * dphi_Ps_CPDs$mean.delta_phi.rad[facet]) 
    dphi_Ps_CPDs$v.mean[facet] <- 1/(sqrt(3) * dphi_Ps_CPDs$mean.delta_phi.rad[facet]) 
    # v = 1/(2 * delta.phi) (square lattice); v= 1/(sqrt(3) * delta.phi) (hexagonal lattice)
    # CPD = 1/(2 * delta.phi) according to Feller Surf and Turf
  }
  
  # # some corrections that are only needed for faulty eyes (so far only damselfly)
  # dphi_Ps_CPDs$CPD.mean[which(dphi_Ps_CPDs$CPD.mean == Inf)] <- mean(dphi_Ps_CPDs$CPD.mean[which(dphi_Ps_CPDs$CPD.mean != Inf)])
  
  
  if(verbose == TRUE){
    cat("Plotting infos to plot device...\n")
  }
  par(mfrow = c(5,1))
  hist(dphi_Ps_CPDs$size, 
       breaks = seq(min(dphi_Ps_CPDs$size, na.rm = TRUE), 
                    max(dphi_Ps_CPDs$size, na.rm = TRUE), 
                    length.out=16),
       main = "Facet diameter",
       xlab = "Facet diameter (um)")
  hist(dphi_Ps_CPDs$mean.delta_phi.deg, 
       breaks = seq(min(dphi_Ps_CPDs$mean.delta_phi.deg, na.rm = TRUE), 
                    max(dphi_Ps_CPDs$mean.delta_phi.deg, na.rm = TRUE), 
                    length.out=16),
       main = "Inter-facet angles",
       xlab = "IF-angle (°)")
  
  hist(dphi_Ps_CPDs$P.mean, 
       breaks = seq(min(dphi_Ps_CPDs$P.mean, na.rm = TRUE), 
                    max(dphi_Ps_CPDs$P.mean, na.rm = TRUE), 
                    length.out=16),
       main = "Eye Paraneter (P)",
       xlab = "P")
  hist(dphi_Ps_CPDs$v.mean, 
       breaks = seq(min(dphi_Ps_CPDs$v.mean, na.rm = TRUE),
                    max(dphi_Ps_CPDs$v.mean, na.rm = TRUE), 
                    length.out=16),
       main = "Acuity (v)",
       xlab = "v")
  hist(dphi_Ps_CPDs$CPD.mean, 
       breaks = seq(min(dphi_Ps_CPDs$CPD.mean, na.rm = TRUE),
                    max(dphi_Ps_CPDs$CPD.mean, na.rm = TRUE), 
                    length.out=16),
       main = "Acuity (CPD)",
       xlab = "CPD")
  par(mfrow = c(1,1))
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("Plotting infos to", plot_file, "\n")
    }
    
    # PDF plots
    pdf(plot_file, # , today()
        onefile = TRUE, paper = "a4")
    
    par(mfrow = c(2,1))
    hist(dphi_Ps_CPDs$size, 
         breaks = seq(min(dphi_Ps_CPDs$size, na.rm = TRUE), 
                      max(dphi_Ps_CPDs$size, na.rm = TRUE), 
                      length.out=16),
         main = "Facet size",
         xlab = "Facet size (um)")
    hist(dphi_Ps_CPDs$mean.delta_phi.deg, 
         breaks = seq(min(dphi_Ps_CPDs$mean.delta_phi.deg, na.rm = TRUE), 
                      max(dphi_Ps_CPDs$mean.delta_phi.deg, na.rm = TRUE), 
                      length.out=16),
         main = "Inter-facet angles",
         xlab = "IF-angle (°)")
    
    par(mfrow = c(3,1))
    hist(dphi_Ps_CPDs$P.mean, 
         breaks = seq(min(dphi_Ps_CPDs$P.mean, na.rm = TRUE), 
                      max(dphi_Ps_CPDs$P.mean, na.rm = TRUE), 
                      length.out=16),
         main = "Eye Paraneter (P)",
         xlab = "P")
    hist(dphi_Ps_CPDs$v.mean, 
         breaks = seq(min(dphi_Ps_CPDs$v.mean, na.rm = TRUE),
                      max(dphi_Ps_CPDs$v.mean, na.rm = TRUE), 
                      length.out=16),
         main = "Acuity (v)",
         xlab = "v")
    hist(dphi_Ps_CPDs$CPD.mean, 
         breaks = seq(min(dphi_Ps_CPDs$CPD.mean, na.rm = TRUE),
                      max(dphi_Ps_CPDs$CPD.mean, na.rm = TRUE), 
                      length.out=16),
         main = "Acuity (CPD)",
         xlab = "CPD")
    par(mfrow = c(1,1))
    
    dev.off()
  }
  
  
  
  return(dphi_Ps_CPDs %>% 
           dplyr::select(facet, mean.delta_phi.deg, mean.delta_phi.rad, P.mean, v.mean, CPD.mean) %>% 
           dplyr::rename(delta_phi.deg = mean.delta_phi.deg,
                         delta_phi.rad = mean.delta_phi.rad,
                         P = P.mean,
                         v = v.mean,
                         CPD = CPD.mean))
}