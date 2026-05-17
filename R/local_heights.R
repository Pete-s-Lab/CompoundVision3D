#' Find Local Vertex Heights
#'
#' Calculate distance of vertices from local plane.
#'
#' @param df A tibble containing triangle center coordinates in columns `x, y, z`.
#' @param search_diam A numerical value defining the size of the search diameter 
#' of defining the local plane.
#' @param cores A numerical value of how many cores to use. Default: `1`.
# @param normalize A numerical value to specifiy if local distances should be
# normalized within a the given radius. If `NULL`, normalizing will be skipped.
# Default: `NULL`. Recommended: `serach_daim/2`.
#'
#' @return Tibble `df` with additional column `local_height`.
#'
#' @export
#' @examples
#' # xxx: add example
#'
calculate_local_heights <- function(df,
                                    search_diam,
                                    cores = 1,
                                    log_scale = TRUE,
                                    plot_file = NULL,
                                    verbose = FALSE,
                                    invert = FALSE){
  
  # load package for multi-core
  require(doParallel)
  
  # # testing
  # # tri_centers_normals <- readr::read_csv("X:/Pub/2021/_Ruehr_AntVision/data/3_triangle_centers_and_normals//AV00001_Camponotus_hyatti_eye1_surface.csv",
  # #                                        show_col_types = FALSE)
  # df = tri_centers_normals[1:100, ]
  # search_diam = curr_facet_estimate*3
  # cores = 12
  # # log_scale = TRUE
  # plot_file = file.path(local_heights_folder,
  #                       gsub("csv$", "pdf", curr_filename_out))
  # verbose = TRUE
  # invert = TRUE
  # #/ testing
  
  # dplyr NULLs
  x <- y <- z <- value <- value.1 <- value.2 <- row_number <- norm.x <- 
    norm.y <- norm.z <- i <- search_diam_local_height <- local_height <- NULL
  
  # convert search_diam to numeric if necessary
  search_diam <- as.numeric(search_diam)
  
  # define function to normalize vector
  normalize_vector <- function(v){
    v/sqrt(sum(v^2))
  }
  
  # calculate once, not once per vertex
  search_diam_local_height <- round(1/8 * search_diam, 2)
  
  # extract coordinates and normals once for faster repeated access
  coords <- as.matrix(df[, c("x", "y", "z")])
  normals <- as.matrix(df[, c("norm.x", "norm.y", "norm.z")])
  
  storage.mode(coords) <- "double"
  storage.mode(normals) <- "double"
  
  # build simple spatial grid to avoid filtering the full data frame for every point
  use_grid <- is.finite(search_diam_local_height) &&
    search_diam_local_height > 0 &&
    all(is.finite(coords))
  
  if(use_grid){
    
    if(verbose == TRUE){
      cat("Building spatial grid...\n")
    }
    
    coord_origin <- apply(coords, 2, min)
    
    cell_mat <- floor(sweep(coords, 2, coord_origin, FUN = "-") /
                        search_diam_local_height)
    
    cell_keys <- paste(cell_mat[, 1],
                       cell_mat[, 2],
                       cell_mat[, 3],
                       sep = "_")
    
    cell_map <- split(seq_len(nrow(df)), cell_keys)
    
    neighbour_offsets <- expand.grid(dx = -1:1,
                                     dy = -1:1,
                                     dz = -1:1,
                                     KEEP.OUT.ATTRS = FALSE)
  }
  
  if(verbose == TRUE){
    cat("Starting analyses on cluster...\n")
    start_time <- Sys.time()
  }
  
  # calculate distance of all vertices to local plane within search_diam
  registerDoParallel(cores)
  
  if(verbose == TRUE){
    cat("Calculating local heights for all", nrow(df), "vertices...\n")
  }
  
  local_heights <- foreach(i = 1:nrow(df),
                           .combine = c,
                           .packages = c('dplyr', 'geometry')) %dopar% {
                             
                             curr.facet.x.y.z <- coords[i, ]
                             
                             if(use_grid){
                               
                               # get all points from the current and directly neighbouring grid cells
                               curr_cell <- cell_mat[i, ]
                               
                               candidate_keys <- paste(curr_cell[1] + neighbour_offsets$dx,
                                                       curr_cell[2] + neighbour_offsets$dy,
                                                       curr_cell[3] + neighbour_offsets$dz,
                                                       sep = "_")
                               
                               idx <- unlist(cell_map[candidate_keys], use.names = FALSE)
                               
                               if(length(idx) > 0){
                                 candidate_coords <- coords[idx, , drop = FALSE]
                                 
                                 # exact same cubic neighbourhood filter as in the original function
                                 keep <- candidate_coords[, 1] >= curr.facet.x.y.z[1] - search_diam_local_height &
                                   candidate_coords[, 2] >= curr.facet.x.y.z[2] - search_diam_local_height &
                                   candidate_coords[, 3] >= curr.facet.x.y.z[3] - search_diam_local_height &
                                   candidate_coords[, 1] <= curr.facet.x.y.z[1] + search_diam_local_height &
                                   candidate_coords[, 2] <= curr.facet.x.y.z[2] + search_diam_local_height &
                                   candidate_coords[, 3] <= curr.facet.x.y.z[3] + search_diam_local_height
                                 
                                 idx <- idx[!is.na(keep) & keep]
                               }
                               
                             } else {
                               
                               # fallback: original full search, but using matrix operations instead of dplyr
                               keep <- coords[, 1] >= curr.facet.x.y.z[1] - search_diam_local_height &
                                 coords[, 2] >= curr.facet.x.y.z[2] - search_diam_local_height &
                                 coords[, 3] >= curr.facet.x.y.z[3] - search_diam_local_height &
                                 coords[, 1] <= curr.facet.x.y.z[1] + search_diam_local_height &
                                 coords[, 2] <= curr.facet.x.y.z[2] + search_diam_local_height &
                                 coords[, 3] <= curr.facet.x.y.z[3] + search_diam_local_height
                               
                               idx <- which(!is.na(keep) & keep)
                             }
                             
                             if(length(idx) == 0){
                               curr_local_height <- NaN
                             } else {
                               
                               # calculate current average facet normal
                               curr.facets.av.normal <- c(mean(normals[idx, 1]), 
                                                          mean(normals[idx, 2]), 
                                                          mean(normals[idx, 3]))
                               
                               # calculate current facet center
                               curr.facets.center <- c(mean(coords[idx, 1]), 
                                                       mean(coords[idx, 2]), 
                                                       mean(coords[idx, 3]))
                               
                               # create unit normal vector of plane normal
                               curr.facets.av.normal.normailzed <- normalize_vector(curr.facets.av.normal)
                               
                               # find vector of current point to arbitrary other point on plane (here: plane center)
                               vector_point_facet_center <- curr.facet.x.y.z - curr.facets.center
                               
                               curr_local_height <- dot(vector_point_facet_center,
                                                        curr.facets.av.normal.normailzed, 
                                                        d = TRUE)
                             }
                             
                             tmp <- curr_local_height
                           }
  
  stopImplicitCluster()
  
  # add distances to local planes to df tibble
  if(invert == FALSE){
    df$local_height <- as.numeric(local_heights)
    # save(df, file = "./data/AV00003_df.Rdata")
  } else{
    df$local_height <- -1 * as.numeric(local_heights)
  }
  
  # save(df, file = "./data/AV00003_df_NORM.Rdata")
  # load(file = "./data/AV00003_df_NORM.Rdata")
  
  
  if(verbose == TRUE){
    cat("Cluster analysis finished.\n")
    end_time <- Sys.time()
    cat("Time difference:", end_time - start_time, "\n")
  }
  
  if(verbose == TRUE){
    cat("Adding quantile-filtered and logarhitmic scales...\n")
  }
  # add color values for RAW local heights
  local_height_cols_raw <- get_height_colors(heights = df$local_height)
  
  
  local_height_cols_filtered <- get_height_colors(heights = df$local_height,
                                                  lower_quantile = 0.1,
                                                  upper_quantile = 0.9)
  
  local_height_cols_filtered_log <- get_height_colors(heights = 10^df$local_height,
                                                      lower_quantile = 0.5,
                                                      upper_quantile = 0.9)
  
  # add colour column for raw values
  df$local_height_col <- local_height_cols_raw
  
  # add colour column for quantile-filtered values
  df$local_height_filterd_col <- local_height_cols_filtered
  
  # add colour and value column for log values
  df$local_height_log <- 10^df$local_height
  df$local_height_log_col <- local_height_cols_filtered_log
  
  
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("2D-plotting to ", plot_file, "...\n")
    }
    # plot height colours for raw, filtered and log-transformed local heights
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    pdf(file = plot_file,
        width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    par(mfrow=c(1,3))
    plot(df$local_height, col = local_height_cols_raw, pch=16, cex=.2,
         main = "raw",
         xlab = "idx",
         ylab = "raw local height")
    
    plot(df$local_height, col = local_height_cols_filtered, pch=16, cex=.2,
         main = "quantile",
         xlab = "idx",
         ylab = "quantile-filtered local height")
    
    plot(df$local_height, col = local_height_cols_filtered_log, pch=16, cex=.2,
         main = "quantile & log",
         xlab = "idx",
         ylab = "quantile-filtered log10(local height)")
    par(mfrow=c(1,1))
    dev.off()
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
  }
  
  if(verbose == TRUE){
    cat("Local heights calculated!\n")
  }
  
  return(df)
}

#' Get colour values !HERE: Change description
#'
#' Calculate distance of vertices from local plane.
#'
#' @param df A tibble containing triangle center coordinates in columns `x, y, z`.
#' @param search_diam A numerical value defining the size of the search diameter 
#' of defining the local plane.
#' @param cores A numerical value of how many cores to use. Default: `1`.
# @param normalize A numerical value to specifiy if local distances should be
# normalized within a the given radius. If `NULL`, normalizing will be skipped.
# Default: `NULL`. Recommended: `serach_daim/2`.
#' @importFrom forceR rescale_to_range
#'
#' @return Tibble `df` with additional column `local_height`.
#'
#' @export
#' @examples
#' # xxx: add example
#'
get_height_colors <- function(heights,
                              lower_quantile = NULL,
                              upper_quantile = NULL,
                              verbose = FALSE){
  
  # # testing
  # lower_quantile = NULL
  # upper_quantile = NULL
  # heights = local_heights$local_height
  # heights = 10^local_heights$local_height
  # range(heights)
  # lower_quantile = 0.1
  # upper_quantile = 0.9
  
  
  if(!is.null(upper_quantile) & !is.null(lower_quantile)){
    # set upper and lower boundary for local heights
    heights[heights >= 
              quantile(heights, upper_quantile)] <-
      quantile(heights, upper_quantile)
    
    if(min(heights) < 0){
      if(verbose == TRUE){
        print("non-logarigthmic")
      }
      heights[heights <= 
                quantile(heights, lower_quantile)] <-
        quantile(heights, lower_quantile)
    } else{
      if(verbose == TRUE){
        print("logarigthmic")
      }
      heights[heights <= 
                quantile(heights, lower_quantile)] <-
        quantile(heights, lower_quantile)
    }
  }
  
  if(verbose == TRUE){
    message("Adding colours for height values...")
  }
  # create color vector
  color_df <- tibble(local_height = round(seq(0, 100000-1, length.out = 100000)), 
                     local_height_col = grey.colors(100000, start=0.0)) %>% 
    distinct(local_height, .keep_all = TRUE)
  
  # add colors
  local_height_cols <- tibble(local_height = heights) %>% 
    mutate(local_height = round(rescale_to_range(local_height, 1, 100000))) %>% 
    left_join(color_df, by = "local_height") %>% 
    pull(local_height_col)
  
  return(local_height_cols)
}


#' Noramlize Local Vertex Heights
#'
#' xxx: update! Calculate distance of vertices from local plane.
#'
#' @param df A tibble containing triangle center coordinates in columns `x, y, z`.
#' @param search_diam A numerical value defining the size of the search diameter 
#' of defining the local plane.
#' @param cores A numerical value of how many cores to use. Default: `1`.
# @param normalize A numerical value to specifiy if local distances should be
# normalized within a the given radius. If `NULL`, normalizing will be skipped.
# Default: `NULL`. Recommended: `serach_daim/2`.
#'
#' @return Tibble `df` with additional column `local_height`.
#'
#' @export
#' @examples
#' # xxx: add example
#'
normalize_local_heights <- function(df,
                                    normalize_diam,
                                    column_to_normalize = "local_height",
                                    cores = 12,
                                    plot_file = NULL,
                                    plot_results = TRUE,
                                    verbose = FALSE){
  
  # # # testing
  # df = local_heights %>% slice(1:1000)
  # column_to_normalize = "local_height_log" # "local_height" "local_height_log"
  # normalize_diam = curr_facet_estimate
  # cores = 12
  # plot_file = file.path(local_heights_normalized_folder,
  #                       gsub("csv$", "pdf", curr_filename_out))
  # verbose = TRUE
  # # df %>%
  # #   select(x, one_of(column_to_normalize)) %>%
  # #   rename(local_height = 2)
  
  
  # convert search_diam to numeric if necessary
  normalize_diam <- as.numeric(normalize_diam)
  
  # The optimized neighborhood search below uses grid cells with
  # cell size = normalize_diam. This requires a positive search diameter.
  if(is.na(normalize_diam) || normalize_diam <= 0){
    stop("normalize_diam must be a positive numeric value.")
  }
  
  # if(!is.null(plot_file)){
  require(rgl)
  # }
  
  require(forceR)
  require(doParallel)
  require(dplyr)
  
  # dplyr NULLs
  x <- y <- z <- value <- value.1 <- value.2 <- row_number <- norm.x <-
    norm.y <- norm.z <- i <- n <- local_height <- local_heights_quantiles_normalized <- NULL
  
  # # plot eye in 'SEM colors'
  # plot3d(df %>%
  #          select(x,y,z),
  #        col = df %>%
  #          pull(local_height_log_col),
  #        aspect = "iso",
  #        size=8)
  
  df <- df %>%
    mutate(n=row_number())
  
  # -------------------------------------------------------------------------
  # OPTIMIZED NORMALIZATION CORE
  # -------------------------------------------------------------------------
  # Original behaviour:
  # For every center point i:
  #   1. find all points inside the x/y/z cube of radius normalize_diam
  #   2. quantile-filter and rescale the selected column within that local cube
  #   3. return one normalized value for every point in that cube
  # Final value per point:
  #   mean of all normalized values assigned to that point from all overlapping cubes
  #
  # The expensive parts in the old version were:
  #   - scanning the complete df once for every point
  #   - building a very large intermediate tibble and then summarising it
  #
  # This version preserves the same logic, but:
  #   - uses a simple 3D grid to get only nearby candidate points
  #   - accumulates sum/count vectors directly
  # -------------------------------------------------------------------------
  
  n_total <- nrow(df)
  
  # limet number of rows for one analsis to get reasonable calulation times
  max_rows = 2500
  starts <- seq(1,n_total,max_rows)
  
  if(verbose == TRUE){
    if(length(starts) > 1){
      cat("Splitting analysis into ", length(starts), " parts.\n")
    }
  }
  
  # Pull frequently used columns into plain vectors.
  # This avoids repeated dplyr::slice(), dplyr::filter(), dplyr::select()
  # inside the inner loop.
  x_vec <- df$x
  y_vec <- df$y
  z_vec <- df$z
  col_to_norm_vec <- df[[column_to_normalize]]
  
  # Build a simple spatial grid.
  # Cell size equals normalize_diam.
  # Therefore, any point within +/- normalize_diam of a center point
  # can only lie in the same or directly adjacent cell along each axis.
  cell_size <- normalize_diam
  
  valid_coord_rows <- is.finite(x_vec) & is.finite(y_vec) & is.finite(z_vec)
  
  if(!any(valid_coord_rows)){
    stop("No rows with finite x/y/z coordinates found.")
  }
  
  x_min <- min(x_vec[valid_coord_rows])
  y_min <- min(y_vec[valid_coord_rows])
  z_min <- min(z_vec[valid_coord_rows])
  
  cell_x <- floor((x_vec - x_min) / cell_size)
  cell_y <- floor((y_vec - y_min) / cell_size)
  cell_z <- floor((z_vec - z_min) / cell_size)
  
  # Create one character key per occupied grid cell.
  # The value stored in each grid cell is the vector of row indices in that cell.
  cell_key <- paste(cell_x, cell_y, cell_z, sep = "_")
  cell_map <- split(seq_len(n_total)[valid_coord_rows], cell_key[valid_coord_rows])
  
  # Offsets for the 27 cells that need to be inspected:
  # center cell plus all adjacent cells in 3D.
  neighbor_offsets <- expand.grid(dx = -1:1,
                                  dy = -1:1,
                                  dz = -1:1)
  neighbor_offsets <- as.matrix(neighbor_offsets)
  
  # Helper: get candidate indices from nearby grid cells,
  # then apply the exact same cube condition as the original dplyr::filter().
  get_cube_neighbor_indices <- function(i){
    
    if(!valid_coord_rows[i]){
      return(integer(0))
    }
    
    curr_cell_x <- cell_x[i]
    curr_cell_y <- cell_y[i]
    curr_cell_z <- cell_z[i]
    
    neighbor_cells <- cbind(curr_cell_x + neighbor_offsets[, "dx"],
                            curr_cell_y + neighbor_offsets[, "dy"],
                            curr_cell_z + neighbor_offsets[, "dz"])
    
    neighbor_keys <- paste(neighbor_cells[, 1],
                           neighbor_cells[, 2],
                           neighbor_cells[, 3],
                           sep = "_")
    
    candidate_idx <- unlist(cell_map[neighbor_keys], use.names = FALSE)
    
    if(length(candidate_idx) == 0){
      return(integer(0))
    }
    
    # Exact cube filtering.
    # This preserves the old condition:
    # x/y/z must each be within +/- normalize_diam of the center point.
    candidate_idx <- candidate_idx[
      x_vec[candidate_idx] >= x_vec[i] - normalize_diam &
        x_vec[candidate_idx] <= x_vec[i] + normalize_diam &
        y_vec[candidate_idx] >= y_vec[i] - normalize_diam &
        y_vec[candidate_idx] <= y_vec[i] + normalize_diam &
        z_vec[candidate_idx] >= z_vec[i] - normalize_diam &
        z_vec[candidate_idx] <= z_vec[i] + normalize_diam
    ]
    
    candidate_idx
  }
  
  # Helper for foreach result combination.
  # Each worker returns only two vectors:
  #   sum_norm   = accumulated normalized values per original row
  #   count_norm = number of contributions per original row
  combine_norm_results <- function(a, b){
    list(sum_norm = a$sum_norm + b$sum_norm,
         count_norm = a$count_norm + b$count_norm)
  }
  
  if(verbose == TRUE){
    cat("Starting analyses on cluster...\n")
    start_time <- Sys.time()
  }
  
  registerDoParallel(cores)
  
  # Instead of one foreach task per point, use one foreach task per block.
  # Inside each block, we loop over the center points and accumulate directly.
  normalized_accumulated <- foreach(k = seq_along(starts),
                                    .combine = combine_norm_results,
                                    .init = list(sum_norm = numeric(n_total),
                                                 count_norm = numeric(n_total)),
                                    .packages = c('forceR')) %dopar% {
                                      
                                      s <- starts[k]
                                      e <- (s + max_rows - 1)
                                      if(e > n_total) e <- n_total
                                      
                                      # Local worker-side accumulators.
                                      # This avoids returning millions of rows.
                                      curr_sum_norm <- numeric(n_total)
                                      curr_count_norm <- numeric(n_total)
                                      
                                      for(i in s:e){
                                        
                                        curr_neighbor_idx <- get_cube_neighbor_indices(i)
                                        
                                        # Same logic as before:
                                        # only process neighborhoods containing more than one point.
                                        if(length(curr_neighbor_idx) > 1){
                                          
                                          curr_values <- col_to_norm_vec[curr_neighbor_idx]
                                          
                                          # set outliers as quantile values
                                          curr_Q <- quantile(curr_values, probs=c(.10, .90), na.rm = FALSE)
                                          curr_Q_min <- curr_Q[1]
                                          curr_Q_max <- curr_Q[2]
                                          
                                          curr_values_quantiles <- pmin(pmax(curr_values, curr_Q_min),
                                                                        curr_Q_max)
                                          
                                          curr_values_normalized <- rescale_to_range(curr_values_quantiles, 0, 1)
                                          
                                          # Old behaviour:
                                          # one normalized value is emitted for every point in the current ROI.
                                          #
                                          # New implementation:
                                          # directly add those emitted values to the target rows.
                                          curr_sum_norm[curr_neighbor_idx] <-
                                            curr_sum_norm[curr_neighbor_idx] + curr_values_normalized
                                          
                                          curr_count_norm[curr_neighbor_idx] <-
                                            curr_count_norm[curr_neighbor_idx] + 1
                                        }
                                      }
                                      
                                      list(sum_norm = curr_sum_norm,
                                           count_norm = curr_count_norm)
                                    }
  
  # terminate cluster
  stopImplicitCluster()
  
  if(verbose == TRUE){
    cat("Multi-threaded analysis finished.\n")
    end_time <- Sys.time()
    cat(end_time - start_time, "\n")
  }
  
  if(verbose == TRUE){
    cat("Summarizing ", sum(normalized_accumulated$count_norm), " results...\n")
  }
  
  # calculate means per facet - i.e. row = n
  # This replaces:
  #
  # df_normalized_raw %>%
  #   group_by(n) %>%
  #   summarise(local_height_norm = mean(local_heights_quantiles_normalized))
  #
  # but produces the same quantity:
  #   local_height_norm = sum of all assigned normalized values /
  #                       number of assigned normalized values
  df_normalized_summarized <- tibble(n = seq_len(n_total),
                                     local_height_norm =
                                       normalized_accumulated$sum_norm /
                                       normalized_accumulated$count_norm)
  
  df_fin <- df %>%
    left_join(df_normalized_summarized, by = "n") %>%
    select(-n)
  
  # -------------------------------------------------------------------------
  # END OPTIMIZED NORMALIZATION CORE
  # -------------------------------------------------------------------------
  
  # fill NA values (points where no other points were int the search radius around them) with median value
  df_fin <- df_fin %>%
    mutate(local_height_norm = ifelse(is.na(local_height_norm), median(local_height_norm, na.rm=TRUE), local_height_norm))
  
  # hist(df_fin$local_height_norm)
  
  if(verbose == TRUE){
    cat("Adding quantile-filtered and logarithmic scales...\n")
  }
  
  local_height_cols_raw_norm <- get_height_colors(heights = df_fin$local_height_norm)
  
  local_height_cols_filtered_norm <- get_height_colors(heights = df_fin$local_height_norm,
                                                       lower_quantile = 0.1,
                                                       upper_quantile = 0.9)
  
  local_height_cols_filtered_log_norm <- get_height_colors(heights = 10^df_fin$local_height_norm)
  
  # add colour column for raw values
  df_fin$local_height_norm_col <- local_height_cols_raw_norm
  
  # add colour and value column for log values
  df_fin$local_height_log_norm <- 10^df_fin$local_height_norm
  df_fin$local_height_log_norm_col <- local_height_cols_filtered_log_norm
  
  
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("2D-plotting to", plot_file, "...\n")
    }
    # plot height colours for raw, filtered and log-transformed local heights
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    pdf(file = plot_file,
        width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    par(mfrow=c(1,3))
    plot(df_fin$local_height, col = df_fin$local_height_col, pch=16, cex=.2,
         main = "raw",
         xlab = "idx",
         ylab = "raw local height")
    
    plot(df_fin$local_height, col = df_fin$local_height_filterd_col, pch=16, cex=.2,
         main = "quantile",
         xlab = "idx",
         ylab = "quantile-filtered local height")
    
    plot(df_fin$local_height, col = df_fin$local_height_log_col, pch=16, cex=.2,
         main = "quantile & log",
         xlab = "idx",
         ylab = "quantile-filtered log10(local height)")
    par(mfrow=c(1,1))
    dev.off()
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    
    plot_file_3D <- gsub("normalized\\.pdf", "normalized.png", plot_file)
    
    if(verbose == TRUE){
      print(paste0("Saving 3D comparison of original and normalized files to ", plot_file_3D, "..."))
    }
    
    # plot eye in 'SEM colors'
    close3d()
    plot3d(df_fin %>%
             select(x,y,z),
           col = df_fin %>%
             pull(local_height_col),
           aspect = "iso",
           size=3)
    
    points3d(df_fin %>%
               select(x,y,z) %>%
               mutate(x = x+max(x)+0.5*diff(range(x))),
             col = df_fin %>%
               pull(local_height_log_col),
             aspect = "iso",
             size=3)
    
    
    points3d(df_fin %>%
               select(x,y,z) %>%
               mutate(z = z-max(z)-0.5*diff(range(z))),
             col = df_fin %>%
               pull(local_height_norm_col),
             aspect = "iso",
             size=3)
    
    points3d(df_fin %>%
               select(x,y,z) %>%
               mutate(x = x+max(x)+0.5*diff(range(x)),
                      z = z-max(z)-0.5*diff(range(z))),
             col = df_fin %>% pull(local_height_log_norm_col),
             aspect = "iso",
             size=3)
    
    # # View along the X-axis
    # view3d(userMatrix = rotationMatrix(pi/2, 0, -1, 0))
    
    par3d(windowRect = c(20, 30, 800, 800))
    # par3d(userMatrix = rotate3d(par3d("userMatrix"), angle2, 0, 0, -1))
    
    # remove bounding box
    # bbox3d(alpha = 0.0, xlab="NULL")
    
    # rotate view roughly to look at eye
    # Step 1: Calculate the mean direction
    mean_vector <- colMeans(df_fin %>%
                              select(norm.x, norm.y, norm.z))
    mean_direction <- mean_vector / sqrt(sum(mean_vector^2))  # Normalize
    
    # Step 2: Calculate the opposite direction
    opposite_direction <- -mean_direction
    
    # # Highlight the mean vector and its opposite
    # arrow3d(rep(0, 3), mean_direction*200, col = "red", lwd = 3)
    # arrow3d(rep(0, 3), opposite_direction*200, col = "green", lwd = 3)
    
    # Step 4: Rotate the view to look in the opposite direction of the mean vector
    view3d(userMatrix = rotationMatrix(acos(-1),  # 180 degrees rotation
                                       opposite_direction[1],
                                       opposite_direction[2],
                                       opposite_direction[3]))
    
    # # rotate view
    # angle1 <- 45 * (pi/180)   # 45 degrees in radians
    # angle2 <- 60 * (pi/180)   # 60 degrees in radians
    # par3d(userMatrix = rotate3d(par3d("userMatrix"), angle1, 0, 1, 0))
    
    rgl.snapshot(plot_file_3D)
    
    if(plot_results != TRUE){
      close3d()
    }
  }
  
  if(verbose == TRUE){
    cat("Normalization done!\n")
  }
  
  return(df_fin)
}



#' Combine facets with landmarks
#'
#' Calculate distance of vertices from local plane.
#'
#' @param df A tibble containing triangle center coordinates in columns `x, y, z`.
#' @param search_diam A numerical value defining the size of the search diameter 
#' of defining the local plane.
#' @param cores A numerical value of how many cores to use. Default: `1`.
# @param normalize A numerical value to specifiy if local distances should be
# normalized within a the given radius. If `NULL`, normalizing will be skipped.
# Default: `NULL`. Recommended: `serach_daim/2`.
#'
#' @return Tibble `df` with additional column `local_height`.
#'
#' @export
#' @examples
#' # xxx: add example
#'
combine_facets_and_LMs <- function(df,
                                   local_heights,
                                   landmarks_file,
                                   crop_log_file,
                                   eyes = 2,
                                   facet_estimate = 14,
                                   cores,
                                   plot_results = FALSE,
                                   verbose = FALSE){
  # # testing
  # df = facet_positions
  # local_heights = local_heights
  # landmarks_file = landmarks_file_name
  # crop_log_file = crop_log_file_name
  # facet_estimate = 14
  # eyes = 4
  # cores = 18
  # plot_results = TRUE
  # verbose = TRUE
  
  # # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  
  # json file loading
  require(rjson)
  
  # load tidyverse for its various conveniences
  require(tidyverse)
  
  # load log file if it exists
  if(!is.na(crop_log_file)){
    if(verbose == TRUE) cat("Loading", crop_log_file, "\n")
    crop_log_data <- read_delim(crop_log_file,
                                delim = " = ",
                                col_names = FALSE,
                                show_col_types = FALSE) %>% 
      dplyr::rename(var = X1,
                    val = X2)
    
    # landmarks ---------------------------------------------------------------
    LMs_df <- read_3DSlicer_landmarks(file = landmarks_file)
  } else{
    warning("There is no *crop.log file in ", crop_log_folder)
    LMs_df <- tibble(LM=character(), x=numeric(),y=numeric(),z=numeric())
    crop_log_data <- tibble(var=character(), val=character())
  }
  
  
  # remove duplicate facets if there are any
  facet_positions <- facet_positions %>% 
    distinct()
  
  
  # local heights -----------------------------------------------------------
  
  # check if units are correct in facet positions file ----------------------
  ranges <- tibble(local_heights = diff(range(local_heights$x)),
                   facet_positions = diff(range(facet_positions$x)),
                   # shrinkwrap = diff(range(shrinkwrap$x)),
                   LMs = diff(range(LMs_df$x)))
  
  if(verbose == TRUE) print(ranges)
  
  flag = 0
  if(ranges$facet_positions[1] <= 0.5 * ranges$local_heights[1]){
    warning("Local height and facets coordinates do not seem to be in the same unit.")
    flag = 1
    
    # correct units if necessary
    facet_positions <- facet_positions %>%
      mutate(x = x*1000,
             y = y*1000,
             z = z*1000)
  }
  
  
  if(ranges$LMs[1] <= 0.5 * ranges$local_heights[1] & curr_CV != "0005"){ 
    # this is true for the three biting midge STLs because the surface was extracted by authors of the existing paper
    # only the Apis mesh (0005) is an exception
    warning("Local height and LM coordinates do not seem to be in the same unit.")
    flag = 1
    
    # correct units if necessary
    LMs_df <- LMs_df %>%
      mutate(x = x*1000,
             y = y*1000,
             z = z*1000)
  }
  
  # check again if units are correct in facet positions file
  ranges <- tibble(local_heights = diff(range(local_heights$x)),
                   facet_positions = diff(range(facet_positions$x)),
                   # shrinkwrap = diff(range(shrinkwrap$x)),
                   LMs = diff(range(LMs_df$x)))
  if(flag == 1){
    if(verbose == TRUE) cat("New ranges:\n")
    if(verbose == TRUE) print(ranges)
  }
  
  
  # plot eye in SEM colors --------------------------------------------------
  # plot eye in 3D to get overview
  if(verbose == TRUE) cat("Plotting 'SEM'-coloured eye of", file_name, "in RGL 3D window.\n")
  if(plot_results == TRUE){
    plot3d(local_heights %>% 
             select(x,y,z), 
           col = local_heights$local_height_log_col, 
           size= 5,
           aspect = "iso")
    title3d(main = "Not yet referenced",
            cex=3)
    
    # plot the facet positions
    spheres3d(facet_positions %>% 
                select(x, y, z), 
              col="red", radius=facet_estimate/4, alpha = 1)
    
    # plot the LM coordinates
    spheres3d(LMs_df %>% 
                select(x, y, z), 
              col="blue", radius=facet_estimate/3, alpha = 1)
    text3d(LMs_df %>% 
             select(x, y, z), 
           texts = LMs_df %>% 
             pull(LM),
           col="blue", radius=facet_estimate/3, alpha = 1)
  }
  
  # create a distance matrix to remove facet candidates that are so close that they shold be treated as duplicates 
  if(verbose == TRUE) cat("Creating distance matrix...\n")
  facet_distance_matrix <- dist(facet_positions %>% 
                                  select(x,y,z),
                                method = "euclidean",
                                diag = TRUE)
  
  # transform distance matrix into data frame
  facet_distance_df_all <- reshape2::melt(as.matrix((facet_distance_matrix), 
                                                    varnames = c("row", "col"))) %>% 
    as_tibble() %>% 
    filter(value > 0)
  
  colnames(facet_distance_df_all) <- c("facet_1", "facet_2", "distance")
  
  
  # check if some facets are very close
  potential_close_facets <- facet_distance_df_all %>% 
    filter(distance <= facet_estimate/3/1.5)
  
  if(nrow(potential_close_facets) > 0){
    if(verbose == TRUE) cat("Removing", nrow(potential_close_facets), "facet candidates removed as duplicates.\n")
    
    # remove facets that may be too close
    close_IDs <- facet_positions %>% 
      slice(potential_close_facets %>% 
              pull(facet_2)) %>% 
      pull(ID) %>% 
      unique() %>% 
      sort()
    
    # filter out facets
    facet_positions_new <- facet_positions %>% 
      filter(ID %in% close_IDs == FALSE)
  } else{
    # no changes
    facet_positions_new <- facet_positions
  }
  
  # plot eye in SEM colors --------------------------------------------------
  # plot eye in 3D to get overview
  if(verbose == TRUE) cat("Plotting 'SEM'-coloured eye of", file_name, "in RGL 3D window.\n")
  if(plot_results == TRUE){
    plot3d(local_heights %>% 
             select(x,y,z), 
           col = local_heights$local_height_log_col, 
           size= 5,
           aspect = "iso")
    
    # plot the facet positions
    spheres3d(facet_positions_new %>% 
                select(x, y, z), 
              col="blue", radius=facet_estimate/4, alpha = 1)
  }
  
  if(nrow(crop_log_data) > 0){
    
    # extract pixel sizes
    # original high resolution scan
    px_size_eyes <- crop_log_data %>% 
      filter(var == "px_size") %>% 
      mutate(val = gsub("(^.+) .+", "\\1", iconv(val, from = "ISO-8859-1", to = "UTF-8"))) %>% 
      pull(val) %>% 
      as.numeric()
    
    # head
    px_size_head <- crop_log_data %>% 
      filter(var == "px_size_head") %>% 
      mutate(val = gsub("(^.+) .+", "\\1", iconv(val, from = "ISO-8859-1", to = "UTF-8"))) %>% 
      pull(val) %>% 
      as.numeric()
    
    # extract x,y coordinates of head crop and eyes 1 and 2 ROIs ---------------------------------------
    
    ROIs <- c("head", paste0("eye", 1:eyes))
    
    
    curr_ROI_coordinates <- get_ROI_coordinates(ROIs,
                                                crop_log_data)
    
    
    # add eye ROI coordinates to curr. eye ROI
    local_heights_translated <- translate_ROIs(df = local_heights,
                                               ROI_coordinates = curr_ROI_coordinates,
                                               eye = as.numeric(curr_eye),
                                               px_size_eyes = px_size_eyes)
    
    facet_positions_translated <- translate_ROIs(df = facet_positions_new,
                                                 ROI_coordinates = curr_ROI_coordinates,
                                                 eye = as.numeric(curr_eye),
                                                 px_size_eyes = px_size_eyes)
  } else{
    warning("No log file found to translate data.")
    local_heights_translated <- local_heights
    facet_positions_translated <- facet_positions
  }
  
  
  # plot eye and LMs in SEM colors to check if trasnlation was successful ------
  if(plot_results == TRUE){
    # plot landmarks
    plot3d(LMs_df %>% 
             select(x, y, z), 
           col="blue", size=10, alpha = 1,
           aspect = "iso")
    
    text3d(LMs_df %>% 
             select(x, y, z),
           texts = LMs_df$LM,
           col="blue")
    
    # plot eye
    points3d(local_heights_translated %>% 
               select(x,y,z), 
             col = local_heights$local_height_log_col, 
             size= 5)
    
    # plot facet positions
    spheres3d(facet_positions_translated %>% 
                select(x, y, z), 
              col="red", radius=facet_estimate/4, alpha = 1)
    
    # combine all to one tibble
    LMs_facets_combined <- rbind(facet_positions_translated %>% 
                                   dplyr::mutate(ID = as.character(ID),
                                                 type = "facet"), 
                                 LMs_df %>% 
                                   dplyr::rename(ID = LM) %>% 
                                   dplyr::mutate(type = "LM"))
    
  }
  if(verbose == TRUE) cat("Facets and landmarks combined!\n")
  
  return(LMs_facets_combined)
}
