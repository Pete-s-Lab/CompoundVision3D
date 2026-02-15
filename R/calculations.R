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
                           eye,
                           px_size_eyes){
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



#' Find facet neighbours and Calculate facet sizes
#'
#' xxx: add description
#'
#' @param df A `tibble` containing facet coordinates in columns `x, y, z`.
#' @param facet_size A `numeric` value containing the estimated facet size.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`. 
#' @return Returns a `tibble` containing the additional columns with info on 
#' facet size, facet neighbours and the number of neighbours of each facet.
#'
#' @export
#' @examples
#' xxx: add example
#'
find_neighbours_deprecated <- function(df,
                                       facet_size,
                                       neighbour_threshold = 1.5,
                                       cores = 1,
                                       plot_results = FALSE,
                                       plot_file = NULL,
                                       verbose = FALSE){
  require(doParallel)
  require(parallel)
  require(doSNOW)
  require(progress)
  require(reshape2)
  
  # # testing
  # df = curr_facets
  # facet_size = facet_estimate
  # neighbour_threshold = 1.5
  # cores = cores
  # plot_results = plot_results
  # plot_file = plot_file
  # verbose = verbose
  
  facet_size <- as.numeric(facet_size)
  
  if(verbose == TRUE) cat("Facet size estimate: ", facet_size, "\n")
  
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
  
  # get original facet IDs
  facet_distance_df$facet_1 <- as.character(df$ID[facet_distance_df$facet_1])
  facet_distance_df$facet_2 <- as.character(df$ID[facet_distance_df$facet_2])
  
  # calculate facet sizes
  if(verbose == TRUE){
    cat(paste0("Calculating temporary sizes for ", nrow(df), " facets according to their closest 2 neighbours (multi-threaded)...\n"))
  }
  
  registerDoParallel(cores)
  facet = "6290"
  df_sizes <- foreach(facet = df %>% 
                        pull(ID) %>% 
                        as.character(), # nrow(df)
                      .combine=rbind, .packages=c('dplyr')) %dopar% {
                        
                        facet_size <- facet_distance_df %>%
                          filter(facet_1 == facet)  %>%
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
  
  # hist(df_sizes)
  
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
  
  facet="6290"
  registerDoParallel(cores)
  neighbour_columns <- foreach(facet = df %>% 
                                 pull(ID) %>% 
                                 as.character(),
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
  
  # df %>% filter(ID == "6290")
  
  # calculate median distance (~facet size) of all facets to their closest three neighbours xxx: this is median of all facets yet - bit maybe enough?
  median.size.tmp <- median(df$size, na.rm = TRUE) * 3
  
  # find neighbours that are closer than median.size.tmp to current facet and remove them from facet.infos
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
    pb$tick(tokens = list(facet = df$ID[n]))
  } 
  
  # df %>% mutate(n = row_number()) %>% filter(ID == "6352")
  
  opts <- list(progress = progress)
  l=5565
  neighbours_tmp <- foreach(l = 1:nrow(df ),# nrow(df)
                            .combine=rbind,
                            .packages=c('dplyr'),
                            .options.snow = opts,
                            .errorhandling = "stop") %dopar% {
                              
                              
                              # registerDoParallel(cores)
                              # l=282
                              # neighbours_tmp <- foreach(l = 1:nrow(df), # nrow(df)
                              #                           .combine=rbind, .packages=c('dplyr')) %dopar% {
                              
                              curr_facet <- df$ID[l]
                              
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
                              
                              # filter actual neighbours of border-facets
                              neighbours_raw_filtered <- neighbours_raw %>% 
                                slice(1:6) %>% 
                                filter(distance <= mean(.$distance[1:2])*1.5)
                              
                              
                              neighbouring_facets <- neighbours_raw_filtered %>%
                                pull(facet_2)
                              
                              # take a max. of 6 and min of 3 neighbours
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
    left_join(neighbours_tmp %>% 
                rename(ID = facet), by="ID")
  
  # df %>% filter(ID == "6290")
  
  if(plot_results == TRUE){
    # test plot
    plot3d(df %>%
             select(x,y,z),
           radius = df$size/2,
           # size = mean(df$size)*1.2,
           col = viridis::viridis(n=6)[df$number.of.neighbours],
           type="s",
           aspect = "iso")
    
    texts3d(df %>%
              select(x,y,z),
            texts = df$ID,
            pos = 1)
    # tmp=1
    # for(tmp in 7:nrow(df)){
    #   plot_facet_and_neighbours(df,
    #                             facet = tmp,
    #                             radius = 2,
    #                             threeD = TRUE)
    #   
    #   invisible(readline(prompt="Press [enter] to continue"))
    # }
  }
  
  # sort(df$number.of.neighbours)
  # hist(df$number.of.neighbours)
  
  # find facet sizes according to their neighbours
  if(verbose == TRUE){
    cat("Finding facet sizes according to closest neighbours (multi-threaded)...\n")
  }
  
  registerDoParallel(cores)
  u=1
  u=352
  eye_sizes <- foreach(u = 1:nrow(df), # nrow(df)
                       .combine=rbind, .packages=c('dplyr', 'filesstrings')) %dopar% {
                         curr_facet <- df$ID[u]
                         
                         curr.neighbours <- as.numeric(str_split(df$neighbours[df$ID == curr_facet], pattern = "; ")[[1]])
                         
                         # if(all(!is.na(curr.neighbours))){
                         curr_mean_distance <- facet_distance_df %>%
                           filter(facet_1 == curr_facet, facet_2 != curr_facet) %>%
                           filter(facet_2 %in% curr.neighbours) %>%
                           mutate(mean_distance = mean(distance)) %>%
                           pull(mean_distance) %>%
                           unique()
                         # } else{
                         #   curr_mean_distance <- 0
                         # }
                         
                         # df$size[df$ID == curr_facet] <- curr_mean_distance
                         tmp <- tibble(ID = curr_facet,
                                       size = curr_mean_distance)
                       }
  stopImplicitCluster()
  
  # hist(eye_sizes$size)
  
  # add sizes to df
  df <- df %>% 
    select(-size) %>% 
    left_join(eye_sizes %>% 
                select(ID, size),
              by="ID")
  
  # cleaning
  df <- df %>% 
    # select(-c(size)) %>% 
    # rename(size = size_final) %>%
    select(ID, neighbours, number.of.neighbours, size)
  
  # calculate mean of sizes according to facet neighbours
  if(verbose == TRUE){
    cat("Calculating mean sizes according to all their actual neighbours...\n")
  }
  
  df$mean_size <- NA
  l=10
  for(l in 1:nrow(df)){
    curr_facet <- df$ID[l]
    curr_neighbours <- as.numeric(str_split(df$neighbours[df$ID == curr_facet], pattern = "; ")[[1]])
    
    # get their sizes
    curr_mean_size <- df %>% 
      filter(ID %in% curr_neighbours) %>% 
      summarise(mean_size = mean(size)) %>% 
      pull(mean_size)
    
    df$mean_size[l] <- curr_mean_size
    
  }
  
  # remove temporary size column
  df <- df %>% 
    select(-size) %>% 
    dplyr::rename(size = mean_size)
  
  if(plot_results == TRUE){
    if(verbose == TRUE){
      cat("Plotting infos to plot device...\n")
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
  }
  
  
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
    cat("Neighbours found!\n")
  }
  
  return(df)
}





#' Find facet neighbours and Calculate facet sizes
#'
#' xxx: add description
#'
#' @param df A `tibble` containing facet coordinates in columns `x, y, z`.
#' @param facet_size A `numeric` value containing the estimated facet size.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @param verbose A `logical` value indicating if message printing is permitted.
#' Default: `FALSE`. 
#' @return Returns a `tibble` containing the additional columns with info on 
#' facet size, facet neighbours and the number of neighbours of each facet.
#'
#' @export
#' @examples
#' xxx: add example
#'
find_neighbours <- function(df,
                            x = "x", y = "y", z = "z",
                            id = "ID",
                            center = TRUE,
                            k_local = 6,
                            knn_search = 20,
                            edge_tol = 0.35,
                            max_neighbours = 6,
                            color_option = "D") {
  
  # # testing
  # df = curr_facets
  # x = "x"
  # y = "y"
  # z = "z"
  # id = "ID"
  # center = TRUE
  # k_local = 6
  # knn_search = 20
  # edge_tol = 0.5
  # max_neighbours = 6
  # color_option = "D"
  
  if (!requireNamespace("geometry", quietly = TRUE))
    stop("Package 'geometry' required. install.packages('geometry')")
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("Package 'RANN' required. install.packages('RANN')")
  if (!requireNamespace("viridisLite", quietly = TRUE))
    stop("Package 'viridisLite' required. install.packages('viridisLite')")
  
  stopifnot(all(c(x, y, z, id) %in% names(df)))
  if (anyDuplicated(df[[id]]) > 0) stop("IDs must be unique.")
  if (knn_search < k_local) stop("knn_search must be >= k_local")
  
  pts <- as.matrix(df[, c(x, y, z)])
  storage.mode(pts) <- "double"
  n <- nrow(pts)
  
  if (center) pts <- sweep(pts, 2, colMeans(pts), "-")
  
  norms <- sqrt(rowSums(pts^2))
  if (any(norms == 0)) stop("Some points have zero norm after centering; cannot normalize.")
  u <- pts / norms
  
  ang_dist <- function(ui, uj) {
    d <- sum(ui * uj)
    d <- max(min(d, 1), -1)
    acos(d)
  }
  
  nn <- RANN::nn2(u, u, k = knn_search + 1)
  idx  <- nn$nn.idx[, -1, drop = FALSE]
  dch  <- nn$nn.dists[, -1, drop = FALSE]
  dang <- 2 * asin(pmin(dch / 2, 1))
  
  local_scale <- apply(dang[, 1:k_local, drop = FALSE], 1, median, na.rm = TRUE)
  
  tri <- geometry::convhulln(u, options = "Qt")
  if (is.null(tri) || nrow(tri) == 0) stop("convhulln failed (degenerate point set?).")
  
  edge_pairs <- function(a, b) cbind(pmin(a, b), pmax(a, b))
  e <- rbind(
    edge_pairs(tri[,1], tri[,2]),
    edge_pairs(tri[,2], tri[,3]),
    edge_pairs(tri[,3], tri[,1])
  )
  e <- unique(e)
  colnames(e) <- c("i", "j")
  
  keep <- logical(nrow(e))
  for (k in seq_len(nrow(e))) {
    i <- e[k, "i"]; j <- e[k, "j"]
    th <- (1 + edge_tol) * max(local_scale[i], local_scale[j])
    keep[k] <- ang_dist(u[i,], u[j,]) <= th
  }
  e <- e[keep, , drop = FALSE]
  
  adj <- vector("list", n)
  if (nrow(e) > 0) {
    for (k in seq_len(nrow(e))) {
      i <- e[k, "i"]; j <- e[k, "j"]
      adj[[i]] <- c(adj[[i]], j)
      adj[[j]] <- c(adj[[j]], i)
    }
    adj <- lapply(adj, unique)
  } else {
    adj <- replicate(n, integer(0), simplify = FALSE)
  }
  
  for (i in seq_len(n)) {
    if (length(adj[[i]]) < min(3, k_local)) {
      cand <- idx[i, ]
      need <- min(3, k_local) - length(adj[[i]])
      add <- setdiff(cand, c(i, adj[[i]]))
      if (length(add) > 0)
        adj[[i]] <- unique(c(adj[[i]], add[seq_len(min(need, length(add)))]))
    }
  }
  
  if (!is.null(max_neighbours)) {
    adj <- lapply(seq_len(n), function(i) {
      nb <- adj[[i]]
      if (length(nb) <= max_neighbours) return(nb)
      dots <- drop(u[nb, , drop = FALSE] %*% u[i, ])
      ord <- order(1 - dots)
      nb[ord[seq_len(max_neighbours)]]
    })
  }
  
  ids <- df[[id]]
  
  # ---- ONLY CHANGE: convert neighbour indices to semicolon string ----
  neighbour_strings <- sapply(adj, function(v) {
    if (length(v) == 0) return("")
    paste(ids[v], collapse = "; ")
  })
  
  df$neighbours   <- neighbour_strings
  df$number.of.neighbours <- lengths(adj)
  
  pal <- viridisLite::viridis(6, option = color_option)
  df$number_of_neighs_cols <- pal[pmin(pmax(df$number.of.neighbours, 1), 6)]
  
  
  
  
  return(df)
}



#' Calculate facet sizes according to their neighbours
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
calculate_facet_size <- function(df,
                                 id_col = "ID",
                                 x_col  = "x",
                                 y_col  = "y",
                                 z_col  = "z",
                                 neighbours_col = "neighbours",
                                 sep = ";",
                                 keep_distances = FALSE) {
  stopifnot(all(c(id_col, x_col, y_col, z_col, neighbours_col) %in% names(df)))
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
  })
  
  coords <- df %>%
    select(
      ID = all_of(id_col),
      x  = all_of(x_col),
      y  = all_of(y_col),
      z  = all_of(z_col)
    ) %>%
    mutate(ID = as.character(ID))
  
  long <- df %>%
    transmute(
      ID = as.character(.data[[id_col]]),
      neighbours_raw = as.character(.data[[neighbours_col]])
    ) %>%
    mutate(neighbours_raw = if_else(is.na(neighbours_raw), "", neighbours_raw)) %>%
    separate_rows(neighbours_raw, sep = fixed(sep)) %>%
    transmute(
      ID,
      neighbour_ID = str_trim(neighbours_raw)
    ) %>%
    filter(neighbour_ID != "")
  
  dists <- long %>%
    left_join(coords, by = "ID") %>%
    left_join(coords, by = c("neighbour_ID" = "ID"), suffix = c("", "_nbr")) %>%
    mutate(
      dist = sqrt((x - x_nbr)^2 + (y - y_nbr)^2 + (z - z_nbr)^2)
    )
  
  out <- dists %>%
    group_by(ID) %>%
    summarise(
      facet_size = mean(dist, na.rm = TRUE),
      n_used = sum(!is.na(dist)),
      .groups = "drop"
    )
  
  if (keep_distances) {
    return(list(summary = out, distances = dists))
  }
  return(out)
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
  # df = df_w_sizes
  # cores = 18
  # plot_file = gsub("_neighbour_and_size_data",
  #                  "_normal_data", 
  #                  plot_file)
  # verbose = TRUE
  
  # get mean coordinate of facets
  eyes_mean <- df %>% 
    mutate(mean_x = mean(x),
           mean_y = mean(y),
           mean_z = mean(z)) %>%
    distinct(mean_x, mean_y, mean_z) 
  # 
  # 
  #   plot3d(df %>%
  #            select(x,y,z),
  #          aspect = "iso")
  #   text3d(df %>%
  #            select(x,y,z),
  #          texts = df %>%
  #            pull(ID))
  #   points3d(eyes_mean, size=20)
  
  
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
      curr_facet <- df$ID[l]
      curr_neighbours <- as.numeric(str_split(df$neighbours[l], pattern = "; ")[[1]])
      
      neighbours_to_keep <- curr_neighbours[which(curr_neighbours %in% df$ID)]
      neighbours_to_remove <- curr_neighbours[which(curr_neighbours %in% df$ID == FALSE)]
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
    cat(paste0("Calculating ", nrow(df), " facet normals according to their neighbours' coordinates (multi-threaded)...\n"))
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
    pb$tick(tokens = list(facet = df$ID[n]))
  } 
  
  opts <- list(progress = progress)
  l=308
  df_normals <- foreach(l = 1:nrow(df),# nrow(df)
                        .combine=rbind, 
                        .packages=c('dplyr', 'stringr', 'CompoundVision3D'), 
                        .options.snow = opts,
                        .errorhandling = "stop") %dopar% {
                          
                          # print(l)
                          curr_facet <- df$ID[l]
                          curr_neighbours <- str_split(df$neighbours[l], pattern = "; ")[[1]]
                          
                          # if(all(!is.na(curr_neighbours))){
                          
                          # get coordinates of current facet and neighbours
                          coords_facet <- df %>% 
                            filter(ID == curr_facet) %>% 
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
                          curr_dists <- tibble(n1 = character(),
                                               n2 = character(),
                                               dist = numeric())
                          n1=1
                          n2=2
                          for(n1 in 1:length(curr_neighbours)){
                            for(n2 in 1:length(curr_neighbours)){
                              curr_neighbour_1 <- curr_neighbours[n1]
                              curr_neighbour_2 <- curr_neighbours[n2]
                              coords_neighbour_1 <- df %>%
                                filter(ID == curr_neighbour_1) %>%
                                select(x,y,z) %>%
                                unlist()
                              coords_neighbour_2 <- df %>%
                                filter(ID == curr_neighbour_2) %>%
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
                            # get neighbours for current triangle
                            curr_neighbour_1 <- curr_dists_clean %>%
                              slice(n) %>%
                              pull(n1)
                            curr_neighbour_2 <- curr_dists_clean %>%
                              slice(n) %>%
                              pull(n2)
                            
                            # print(paste0("Building triangle with facets ", curr_facet, ", ", curr_neighbour_1, " and ", curr_neighbour_2))
                            coords_neighbour_1 <- df %>%
                              filter(ID == curr_neighbour_1) %>%
                              select(x,y,z) %>%
                              unlist()
                            coords_neighbour_2 <- df %>%
                              filter(ID == curr_neighbour_2) %>%
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
                          
                          tmp <- tibble(ID = curr_facet,
                                        norm.x = mean(curr_normals_x[1]), 
                                        norm.y = mean(curr_normals_y[1]), 
                                        norm.z = mean(curr_normals_z[1]))
                          # print("******************************")
                        }
  
  stopCluster(cl) 
  
  
  # average facet angles
  l=1
  df_normals_avg <- tibble(ID = character(),
                           norm.x = numeric(),
                           norm.y = numeric(),
                           norm.z = numeric())
  
  for(l in 1:nrow(df)) {
    
    # print(l)
    curr_facet <- df$ID[l]
    curr_neighbours <- str_split(df$neighbours[l], pattern = "; ")[[1]]
    
    curr_normal <- df_normals %>% 
      filter(ID == curr_facet) %>% 
      select(norm.x,norm.y,norm.z)
    
    curr_neighbor_normals <- df_normals %>% 
      filter(ID %in% curr_neighbours) %>% 
      select(norm.x,norm.y,norm.z)
    
    curr_normal_avg_x = mean(c(rep(curr_normal$norm.x, 1),
                               curr_neighbor_normals$norm.x))
    curr_normal_avg_y = mean(c(rep(curr_normal$norm.y, 1),
                               curr_neighbor_normals$norm.y))
    curr_normal_avg_z = mean(c(rep(curr_normal$norm.z, 1),
                               curr_neighbor_normals$norm.z))
    
    df_normals_avg <- df_normals_avg %>% 
      add_row(ID = curr_facet,
              norm.x = curr_normal_avg_x,
              norm.y = curr_normal_avg_y,
              norm.z = curr_normal_avg_z)
  }
  
  # replace original values
  df_normals <- df_normals_avg
  
  if(plot_results == TRUE){
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
  }
  
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
    cat("Normals calculated!\n")
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
                                 plot_results = FALSE,
                                 plot_file = NULL,
                                 verbose = FALSE){
  
  require(doParallel)
  
  # # testing
  # df = df_w_normals
  # cores = 18
  # plot_results = plot_results
  # plot_file = gsub("_neighbour_and_size_data",
  #                  "_optics_parameters",
  #                  plot_file)
  # verbose = TRUE
  
  if(verbose == TRUE){
    cat(paste0("Calculating IF angles, P, v, and CPD for ", nrow(df), " facets (multi-threaded)...\n"))
  }
  
  
  registerDoParallel(cores)
  l=7188
  dphi_Ps_CPDs <- foreach(l = 1:nrow(df), # nrow(df)
                          .combine=rbind, .packages=c('dplyr', 'filesstrings', 'CompoundVision3D')) %dopar% {
                            
                            facet_no <- df$ID[l]
                            
                            curr_facet_coords <- df %>%
                              filter(ID == facet_no) %>%
                              select(x, y, z)
                            
                            curr_facet_normal <- df %>%
                              filter(ID == facet_no) %>%
                              select(norm.x, norm.y, norm.z)
                            
                            if(all(!is.na(curr_facet_normal))){
                              if(all(curr_facet_normal != 0)){
                                # get neighbouring facets without NAs
                                curr_facet_neighbours <- str_split(df %>%
                                                                     filter(ID==facet_no) %>%
                                                                     pull(neighbours),
                                                                   pattern = "; ")[[1]]
                                
                                if(all(!is.na(curr_facet_neighbours))){
                                  
                                  delta_phis.rad <- c()
                                  delta_phis.deg <- c()
                                  
                                  neighbour <- curr_facet_neighbours[1]
                                  for (neighbour in curr_facet_neighbours) {
                                    curr_neighbour_coords <- df %>%
                                      filter(ID == neighbour) %>%
                                      select(x, y, z)
                                    
                                    curr_neighbour_normal <- df %>%
                                      filter(ID == neighbour) %>%
                                      select(norm.x, norm.y, norm.z)
                                    
                                    
                                    curr_delta_phi.rad <- calc_delta.phi(curr_facet_normal, 
                                                                         curr_neighbour_normal, 
                                                                         type = "r")
                                    if(is.na(curr_delta_phi.rad)){
                                      curr_delta_phi.rad = data.frame(0)
                                    }
                                    curr_delta_phi.deg <- curr_delta_phi.rad*180/pi
                                    
                                    
                                    delta_phis.rad <- c(delta_phis.rad, curr_delta_phi.rad[1,1])
                                    delta_phis.deg <- c(delta_phis.deg, curr_delta_phi.deg[1,1])
                                    
                                    # print(paste0("curr. ° = ", curr_delta_phi.deg))
                                  }
                                  
                                  # df$delta_phi.rad[l] <- mean(delta_phis.rad)
                                  # df$delta_phi.deg[l] <- mean(delta_phis.deg)
                                  
                                  curr_P <- df$size[df$ID == as.character(facet_no)] * 
                                    mean(delta_phis.rad)  * (sqrt(3)/2) # eye parameter Snyder 1977. Brigitte: (sqrt(3)/2) *
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
                            
                            tmp <- rbind(tibble(ID = facet_no, 
                                                delta_phi.rad =  mean(delta_phis.rad), 
                                                delta_phi.deg = mean(delta_phis.deg), 
                                                P = curr_P, CPD = curr_CPD))
                            # df$P[l] <- curr_P
                            # df$CPD[l] <- curr_CPD
                          }
  stopImplicitCluster()
  
  # add optic parameters to df
  dphi_Ps_CPDs <- df %>% 
    left_join(dphi_Ps_CPDs, 
              by = "ID")
  
  if(verbose == TRUE){
    cat("Calculating weighted mean IF angles...\n")
  }
  
  # find mean IF angle for each facet
  dphi_Ps_CPDs$mean.delta_phi.deg <- NA
  q = 1
  for(q in 1:nrow(dphi_Ps_CPDs)){
    curr_facet <- dphi_Ps_CPDs$ID[q]
    curr.neighbours <- as.numeric(
      str_split(dphi_Ps_CPDs %>% 
                  filter(ID==curr_facet) %>% 
                  pull(neighbours), 
                pattern = "; ")[[1]])
    
    # get number of neghbors for weighted averaging
    no_of_neihbors <- length(curr.neighbours)
    
    weighting_factor <- 6 - no_of_neihbors
    
    # get mean of IF-angles with curr facet angle doubled-weighted
    dphi_Ps_CPDs$mean.delta_phi.deg[q] <- mean(c(rep(dphi_Ps_CPDs$delta_phi.deg[dphi_Ps_CPDs$ID == curr_facet], (weighting_factor+1)), 
                                                 dphi_Ps_CPDs$delta_phi.deg[dphi_Ps_CPDs$ID == curr_facet],
                                                 dphi_Ps_CPDs$delta_phi.deg[dphi_Ps_CPDs$ID %in% curr.neighbours]))
  }
  dphi_Ps_CPDs$mean.delta_phi.rad <- dphi_Ps_CPDs$mean.delta_phi.deg*pi/180
  
  # calculate P and CPD with average angles per facet
  if(verbose == TRUE){
    cat("Calculating mean P and CPD angles...\n")
  }
  dphi_Ps_CPDs$P.mean <- NA
  dphi_Ps_CPDs$CPD.mean <- NA
  dphi_Ps_CPDs$v.mean <- NA
  f=1
  for(f in 1:nrow(dphi_Ps_CPDs)){
    facet = dphi_Ps_CPDs$ID[f]
    dphi_Ps_CPDs$P.mean[dphi_Ps_CPDs$ID == facet] <- dphi_Ps_CPDs$size[dphi_Ps_CPDs$ID == facet] * 
      dphi_Ps_CPDs$mean.delta_phi.rad[dphi_Ps_CPDs$ID == facet] * 
      (sqrt(3)/2) # eye parameter Snyder 1977. Brigitte: (sqrt(3)/2) * 
    dphi_Ps_CPDs$CPD.mean[dphi_Ps_CPDs$ID == facet] <- 1 / (2 * dphi_Ps_CPDs$mean.delta_phi.rad[dphi_Ps_CPDs$ID == facet]) 
    dphi_Ps_CPDs$v.mean[dphi_Ps_CPDs$ID == facet] <- 1/(sqrt(3) * dphi_Ps_CPDs$mean.delta_phi.rad[dphi_Ps_CPDs$ID == facet]) 
    # v = 1/(2 * delta.phi) (square lattice); v= 1/(sqrt(3) * delta.phi) (hexagonal lattice)
    # CPD = 1/(2 * delta.phi) according to Feller Surf and Turf
  }
  
  # # some corrections that are only needed for faulty eyes (so far only damselfly)
  # dphi_Ps_CPDs$CPD.mean[which(dphi_Ps_CPDs$CPD.mean == Inf)] <- mean(dphi_Ps_CPDs$CPD.mean[which(dphi_Ps_CPDs$CPD.mean != Inf)])
  
  if(plot_results == TRUE){
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
  }
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
  
  
  
  if(verbose == TRUE){
    cat("Optic parameters calcualted!\n")
  }
  
  return(dphi_Ps_CPDs %>% 
           dplyr::select(ID, mean.delta_phi.deg, mean.delta_phi.rad, P.mean, v.mean, CPD.mean) %>% 
           dplyr::rename(delta_phi.deg = mean.delta_phi.deg,
                         delta_phi.rad = mean.delta_phi.rad,
                         P = P.mean,
                         v = v.mean,
                         CPD = CPD.mean))
}






#' Calculates the area of a hexagon given the long diagonal
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
calculate_hexagon_area <- function(long_diagonal) {
  # Calculate the side length
  side_length <- long_diagonal / 2
  
  # Calculate the area using the formula
  area <- (3 * sqrt(3) / 2) * (side_length^2)
  
  return(area)
}




#' Converts 3D coordinates to latitude and longitude
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
convert_to_latlon <- function(x, y, z) {
  # # testing
  # x = curr_df_all %>%
  #   pull(corn.proj.x)
  # y = curr_df_all %>%
  #   pull(corn.proj.y)
  # z = curr_df_all %>%
  #   pull(corn.proj.z)
  
  # Calculate the radius (should be 1 for unit sphere)
  r <- sqrt(x^2 + y^2 + z^2)
  
  # Calculate latitude in radians
  latitude <- asin(z / r)
  
  # Calculate longitude in radians
  longitude <- atan2(y, x)
  
  # Convert latitude and longitude to degrees
  latitude <- latitude * 180 / pi
  longitude <- longitude * 180 / pi
  
  lat_lon <- tibble(latitude,
                    longitude)
  
  return(lat_lon)
}


#' Calculates intersection of vector and sphere = corneal projection
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
vector.sphere.intersect <- function(point, vector, sphere.c, sphere.r){
  # expected formats:
  # point: vector of length 3 with x, y, and z coordinates of starting point (lens)
  # vector: vector of length 3 with vector direction in x, y, and z
  # sphere.c: vector of length 3 with x, y, and z coordinates of center of corneal projection sphere (e.g., c(0, 2571, 975))
  # sphere.r: single number describing radius of corneal projection sphere (e.g. 7750)
  
  p1.x = point[1]
  p1.y = point[2]
  p1.z = point[3]
  
  v1.x = vector[1]
  v1.y = vector[2]
  v1.z = vector[3]
  
  p2.x = p1.x + v1.x
  p2.y = p1.y + v1.y
  p2.z = p1.z + v1.z
  
  sphere.c.x = sphere.c[1]
  sphere.c.y = sphere.c[2]
  sphere.c.z =sphere.c[3]
  
  A = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
  B = 2.0 * (p1.x * v1.x + p1.y * v1.y + p1.z * v1.z - v1.x * sphere.c.x - v1.y * sphere.c.y - v1.z * sphere.c.z);
  C = p1.x * p1.x - 2 * p1.x * sphere.c.x + sphere.c.x * sphere.c.x + p1.y * p1.y - 2 * p1.y * sphere.c.y + sphere.c.y * sphere.c.y +
    p1.z * p1.z - 2 * p1.z * sphere.c.z + sphere.c.z * sphere.c.z - sphere.r * sphere.r;
  
  # discriminant
  D = B * B - 4 * A * C;
  
  t1 = ( -B - sqrt ( D ) ) / ( 2.0 * A );
  
  solution1 = c( p1.x * ( 1 - t1 ) + t1 * p2.x,
                 p1.y * ( 1 - t1 ) + t1 * p2.y,
                 p1.z * ( 1 - t1 ) + t1 * p2.z )
  
  t2 = ( -B + sqrt ( D ) ) / ( 2.0 * A )
  solution2 = c( p1.x * ( 1 - t2 ) + t2 * p2.x,
                 p1.y * ( 1 - t2 ) + t2 * p2.y,
                 p1.z * ( 1 - t2 ) + t2 * p2.z )
  
  return(as_tibble(rbind(solution1, solution2)))
}