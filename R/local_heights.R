#' Calculate Local Surface Heights
#'
#' Calculates the local height of each triangle centre in a surface mesh
#' relative to its surrounding surface. For each triangle, a local reference
#' plane is defined from the mean centre coordinates and mean surface normal of
#' neighbouring triangles. The signed distance of the focal triangle centre
#' from this plane is returned as its local surface height.
#'
#' Local surface heights are used by CV3D to identify the raised centres of
#' compound-eye facets.
#'
#' @param df A data frame or tibble containing triangle-centre coordinates
#'   (`x`, `y`, `z`) and triangle-normal components (`norm.x`, `norm.y`,
#'   `norm.z`).
#' @param neighbourhood_radius Numeric. Radius, in micrometres (µm), of the
#'   spherical Euclidean neighbourhood used to define the local reference plane.
#'   CV3D currently assumes all input coordinates are expressed in micrometres.
#' @param colour_lower_quantile Numeric. Lower quantile used to clip the
#'   greyscale colour representation of raw local heights. Default: `0.10`.
#' @param colour_upper_quantile Numeric. Upper quantile used to clip the
#'   greyscale colour representation of raw local heights. Default: `0.90`.
#' @param contrast_lower_quantile Numeric. Lower quantile used after the
#'   historical `10^local_height` peak-enhancing transform before rescaling the
#'   contrast variable to 0--1. Default: `0.50`.
#' @param contrast_upper_quantile Numeric. Upper quantile used after the
#'   peak-enhancing transform before rescaling to 0--1. Default: `0.90`.
#' @param cores Integer. Number of processor cores used for parallel
#'   calculation. Default: `1`.
#' @param plot_file Character or `NULL`. If a file path is supplied, write a
#'   PDF containing diagnostic plots of the raw, quantile-filtered, and
#'   contrast-enhanced local-height values. Default: `NULL`.
#' @param verbose Logical. If `TRUE`, print progress and timing information.
#'   Default: `FALSE`.
#' @param invert Logical. If `TRUE`, multiply calculated local heights by
#'   `-1`. Default: `FALSE`.
#'
#' @return The input data with five additional columns: `local_height`, the
#'   signed local surface height in micrometres; `local_height_col`, colours representing the
#'   raw local heights; `local_height_filtered_col`, colours based on
#'   quantile-filtered local heights; `local_height_contrast`, a 0--1
#'   peak-enhancing contrast scale based on `10^local_height` after clipping
#'   to configurable quantiles (50th and 90th percentiles by default); and `local_height_contrast_col`, colours
#'   based on that contrast scale.
#'
#' @export
#' @examples
#' surface <- expand.grid(x = -1:1, y = -1:1)
#' surface$z <- 0.1 * (surface$x^2 + surface$y^2)
#' surface$norm.x <- 0
#' surface$norm.y <- 0
#' surface$norm.z <- 1
#' heights <- calculate_local_heights(
#'   surface,
#'   neighbourhood_radius = 1,
#'   cores = 1,
#'   verbose = FALSE
#' )
#' heights[, c("x", "y", "z", "local_height")]
#'
calculate_local_heights <- function(df,
                                    neighbourhood_radius,
                                    colour_lower_quantile = 0.10,
                                    colour_upper_quantile = 0.90,
                                    contrast_lower_quantile = 0.50,
                                    contrast_upper_quantile = 0.90,
                                    cores = 1,
                                    plot_file = NULL,
                                    verbose = FALSE,
                                    invert = FALSE){
  
  # load package for multi-core
  
  # dplyr NULLs
  x <- y <- z <- value <- value.1 <- value.2 <- row_number <- norm.x <- 
    norm.y <- norm.z <- i <- local_height <- NULL
  
  # convert neighbourhood radius to numeric if necessary
  neighbourhood_radius <- as.numeric(neighbourhood_radius)
  if(length(neighbourhood_radius) != 1 || !is.finite(neighbourhood_radius) || neighbourhood_radius <= 0){
    stop("neighbourhood_radius must be a positive finite numeric value.", call. = FALSE)
  }

  validate_quantile_pair <- function(lower, upper, label) {
    if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1 ||
        length(upper) != 1 || !is.finite(lower) || !is.finite(upper) ||
        lower < 0 || upper > 1 || lower >= upper) {
      stop(label, " must satisfy 0 <= lower < upper <= 1.", call. = FALSE)
    }
  }
  validate_quantile_pair(colour_lower_quantile, colour_upper_quantile,
                         "Colour quantiles")
  validate_quantile_pair(contrast_lower_quantile, contrast_upper_quantile,
                         "Contrast quantiles")
  
  # define function to normalize vector
  normalize_vector <- function(v){
    v/sqrt(sum(v^2))
  }
  
  
  # extract coordinates and normals once for faster repeated access
  coords <- as.matrix(df[, c("x", "y", "z")])
  normals <- as.matrix(df[, c("norm.x", "norm.y", "norm.z")])
  
  storage.mode(coords) <- "double"
  storage.mode(normals) <- "double"
  
  # build simple spatial grid to avoid filtering the full data frame for every point
  use_grid <- is.finite(neighbourhood_radius) &&
    neighbourhood_radius > 0 &&
    all(is.finite(coords))
  
  if(use_grid){
    
    if(verbose == TRUE){
      cat("Building spatial grid...\n")
    }
    
    coord_origin <- apply(coords, 2, min)
    
    cell_mat <- floor(sweep(coords, 2, coord_origin, FUN = "-") /
                        neighbourhood_radius)
    
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
  
  # calculate distance of all vertices to local plane within the local neighbourhood
  doParallel::registerDoParallel(cores)
  
  if(verbose == TRUE){
    cat("Calculating local heights for all", nrow(df), "vertices...\n")
  }
  
  local_heights <- foreach::foreach(i = 1:nrow(df),
                           .combine = c) %dopar% {
                             
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
                                 
                                 # Keep only points inside the spherical Euclidean neighbourhood.
                                 # Squared distances avoid an unnecessary square-root operation.
                                 offsets <- sweep(candidate_coords, 2, curr.facet.x.y.z, FUN = "-")
                                 keep <- rowSums(offsets^2) <= neighbourhood_radius^2
                                 
                                 idx <- idx[!is.na(keep) & keep]
                               }
                               
                             } else {
                               
                               # Fallback full search using the same spherical Euclidean criterion.
                               offsets <- sweep(coords, 2, curr.facet.x.y.z, FUN = "-")
                               keep <- rowSums(offsets^2) <= neighbourhood_radius^2
                               
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
                               
                               curr_local_height <- sum(vector_point_facet_center * curr.facets.av.normal.normailzed)
                             }
                             
                             tmp <- curr_local_height
                           }
  
  doParallel::stopImplicitCluster()
  
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
    cat("Adding quantile-filtered and contrast-enhanced scales...\n")
  }
  # add colour values for raw local heights
  local_height_cols_raw <- get_height_colors(heights = df$local_height)

  local_height_cols_filtered <- get_height_colors(
    heights = df$local_height,
    lower_quantile = colour_lower_quantile,
    upper_quantile = colour_upper_quantile
  )

  # Peak-enhancing contrast scale used for visual inspection and thresholding.
  # This retains the historical 10^height transformation and configurable
  # quantile clipping, but exposes the result on an intuitive 0--1 scale.
  exp10_height <- 10^df$local_height
  finite_exp10 <- is.finite(exp10_height)
  local_height_contrast <- rep(NA_real_, length(exp10_height))
  if (any(finite_exp10)) {
    contrast_bounds <- stats::quantile(exp10_height[finite_exp10],
                                       probs = c(contrast_lower_quantile, contrast_upper_quantile),
                                       na.rm = TRUE,
                                       names = FALSE)
    contrast_clipped <- pmin(pmax(exp10_height[finite_exp10], contrast_bounds[1]),
                             contrast_bounds[2])
    if (isTRUE(all.equal(contrast_bounds[1], contrast_bounds[2]))) {
      local_height_contrast[finite_exp10] <- 0.5
    } else {
      local_height_contrast[finite_exp10] <-
        (contrast_clipped - contrast_bounds[1]) / diff(contrast_bounds)
    }
  }
  local_height_cols_contrast <- get_height_colors(heights = local_height_contrast)

  df$local_height_col <- local_height_cols_raw
  df$local_height_filtered_col <- local_height_cols_filtered
  df$local_height_contrast <- local_height_contrast
  df$local_height_contrast_col <- local_height_cols_contrast
  
  
  
  if(!is.null(plot_file)){
    if(verbose == TRUE){
      cat("2D-plotting to ", plot_file, "...\n")
    }
    # plot raw, clipped-colour, and contrast-enhanced local heights
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    grDevices::pdf(file = plot_file,
        width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    graphics::par(mfrow=c(1,3))
    graphics::plot(df$local_height, col = local_height_cols_raw, pch=16, cex=.2,
         main = "raw",
         xlab = "idx",
         ylab = "raw local height")
    
    graphics::plot(df$local_height, col = local_height_cols_filtered, pch=16, cex=.2,
         main = "quantile",
         xlab = "idx",
         ylab = "quantile-filtered local height")
    
    graphics::plot(df$local_height_contrast, col = local_height_cols_contrast, pch=16, cex=.2,
         main = "contrast",
         xlab = "idx",
         ylab = "local-height contrast (0-1)")
    graphics::par(mfrow=c(1,1))
    grDevices::dev.off()
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
  }
  
  if(verbose == TRUE){
    cat("Local heights calculated!\n")
  }
  
  return(df)
}

#' Map Surface Heights to Greyscale Colours
#'
#' Maps numeric surface-height values to a greyscale colour vector for
#' visualisation. Optionally, values can be clipped at lower and upper
#' quantiles before colour mapping to reduce the influence of extreme values.
#'
#' @param heights Numeric vector of surface-height values.
#' @param lower_quantile Numeric value between `0` and `1`, or `NULL`. If
#'   supplied together with `upper_quantile`, values below this quantile are
#'   clipped to the corresponding quantile value. Default: `NULL`.
#' @param upper_quantile Numeric value between `0` and `1`, or `NULL`. If
#'   supplied together with `lower_quantile`, values above this quantile are
#'   clipped to the corresponding quantile value. Default: `NULL`.
#' @param verbose Logical. If `TRUE`, print progress information. Default:
#'   `FALSE`.
#'
#' @return A character vector of greyscale colours with the same length as
#'   `heights`.
#'
#' @keywords internal
#'
get_height_colors <- function(heights,
                              lower_quantile = NULL,
                              upper_quantile = NULL,
                              verbose = FALSE) {
  heights <- as.numeric(heights)

  if (!is.null(lower_quantile) || !is.null(upper_quantile)) {
    if (is.null(lower_quantile) || is.null(upper_quantile)) {
      stop("lower_quantile and upper_quantile must either both be supplied or both be NULL.", call. = FALSE)
    }
    if (!is.numeric(lower_quantile) || !is.numeric(upper_quantile) ||
        length(lower_quantile) != 1 || length(upper_quantile) != 1 ||
        !is.finite(lower_quantile) || !is.finite(upper_quantile) ||
        lower_quantile < 0 || upper_quantile > 1 || lower_quantile > upper_quantile) {
      stop("Quantiles must satisfy 0 <= lower_quantile <= upper_quantile <= 1.", call. = FALSE)
    }
    bounds <- stats::quantile(heights, probs = c(lower_quantile, upper_quantile), na.rm = TRUE)
    heights <- pmin(pmax(heights, bounds[1]), bounds[2])
  }

  if (verbose) message("Adding colours for height values...")

  finite <- is.finite(heights)
  colours <- rep(NA_character_, length(heights))
  if (!any(finite)) return(colours)

  values <- heights[finite]
  value_range <- range(values)
  if (diff(value_range) == 0) {
    scaled <- rep(50000L, length(values))
  } else {
    scaled <- round(1 + (values - value_range[1]) / diff(value_range) * 99999)
  }
  palette <- grDevices::grey.colors(100000, start = 0.0)
  colours[finite] <- palette[pmax(1L, pmin(100000L, as.integer(scaled)))]
  colours
}


#' Normalize Local Surface Heights
#'
#' Normalizes local surface-height values relative to their spatial
#' neighbourhood. For each triangle centre, values within a spherical
#' Euclidean neighbourhood are clipped to configurable lower and upper
#' quantiles (10th and 90th percentiles by default) and rescaled to the range
#' from 0 to 1. Because neighbourhoods overlap, each
#' triangle can receive normalized values from multiple neighbourhoods. The
#' final normalized value of each triangle is the mean of these values.
#'
#' This optional normalization can reduce spatial variation in local-height
#' magnitude across a compound-eye surface before facet-candidate detection.
#'
#' @param df A data frame or tibble containing triangle-centre coordinates
#'   (`x`, `y`, `z`) and the numeric column to be normalized.
#' @param neighbourhood_radius Numeric. Radius, in micrometres (µm), of the
#'   spherical Euclidean neighbourhood used for local normalization. CV3D
#'   currently assumes all input coordinates are expressed in micrometres.
#' @param lower_quantile Numeric. Lower local quantile used for clipping before
#'   0--1 rescaling. Default: `0.10`.
#' @param upper_quantile Numeric. Upper local quantile used for clipping before
#'   0--1 rescaling. Default: `0.90`.
#' @param column_to_normalize Character. Name of the numeric column to
#'   normalize. Default: `"local_height"`.
#' @param cores Integer. Number of processor cores used for parallel
#'   calculation. Default: `1`.
#' @param plot_file Character or `NULL`. If supplied, write diagnostic
#'   visualizations of the original and normalized local-height values.
#'   Default: `NULL`.
#' @param plot_results Logical. If `TRUE` and `plot_file` is supplied, keep the
#'   interactive `rgl` visualization open after saving the diagnostic image.
#'   Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress and timing information.
#'   Default: `FALSE`.
#'
#' @return The input data with additional normalization columns, including
#'   `local_height_norm`, containing locally normalized values between 0 and 1;
#'   `local_height_norm_col`, its greyscale colour representation;
#'   `local_height_norm_contrast`, a base-10 peak-enhancing transform rescaled
#'   to 0--1; and `local_height_norm_contrast_col`, its colour representation.
#'
#' @examples
#' surface <- expand.grid(x = -1:1, y = -1:1)
#' surface$z <- 0
#' surface$local_height <- seq_len(nrow(surface))
#' normalized <- normalize_local_heights(
#'   surface,
#'   neighbourhood_radius = 1.5,
#'   cores = 1,
#'   plot_results = FALSE,
#'   verbose = FALSE
#' )
#' normalized[, c("x", "y", "z", "local_height_norm")]
#'
#' @export
#'
normalize_local_heights <- function(df,
                                    neighbourhood_radius,
                                    column_to_normalize = "local_height",
                                    lower_quantile = 0.10,
                                    upper_quantile = 0.90,
                                    cores = 1,
                                    plot_file = NULL,
                                    plot_results = FALSE,
                                    verbose = FALSE){
  
  # # # testing
  # df = local_heights %>% dplyr::slice(1:1000)
  # neighbourhood_radius = curr_facet_estimate
  # cores = 12
  # plot_file = file.path(local_heights_normalized_folder,
  #                       gsub("csv$", "pdf", curr_filename_out))
  # verbose = TRUE
  # # df %>%
  # #   dplyr::select(x, one_of(column_to_normalize)) %>%
  # #   rename(local_height = 2)
  
  
  # convert neighbourhood radius to numeric if necessary
  neighbourhood_radius <- as.numeric(neighbourhood_radius)
  
  # The optimized neighbourhood search below uses grid cells with
  # cell size = neighbourhood_radius. This requires a positive radius.
  if(length(neighbourhood_radius) != 1 || !is.finite(neighbourhood_radius) || neighbourhood_radius <= 0){
    stop("neighbourhood_radius must be a positive finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(lower_quantile) || !is.numeric(upper_quantile) ||
      length(lower_quantile) != 1 || length(upper_quantile) != 1 ||
      !is.finite(lower_quantile) || !is.finite(upper_quantile) ||
      lower_quantile < 0 || upper_quantile > 1 || lower_quantile >= upper_quantile) {
    stop("Normalization quantiles must satisfy 0 <= lower_quantile < upper_quantile <= 1.", call. = FALSE)
  }
  # dplyr NULLs
  x <- y <- z <- value <- value.1 <- value.2 <- row_number <- norm.x <-
    norm.y <- norm.z <- i <- n <- local_height <- local_heights_quantiles_normalized <- NULL
  
  # # plot eye in 'SEM colors'
  # rgl::plot3d(df %>%
  #          dplyr::select(x,y,z),
  #        col = df %>%
  #        aspect = "iso",
  #        size=8)
  
  df <- df %>%
    dplyr::mutate(n=dplyr::row_number())
  
  # -------------------------------------------------------------------------
  # OPTIMIZED NORMALIZATION CORE
  # -------------------------------------------------------------------------
  # Original behaviour:
  # For every center point i:
  #   1. find all points inside the spherical Euclidean neighbourhood
  #   2. quantile-filter and rescale the selected column within that local sphere
  #   3. return one normalized value for every point in that sphere
  # Final value per point:
  #   mean of all normalized values assigned to that point from all overlapping spheres
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

  rescale_local <- function(values, lower = 0, upper = 1) {
    values <- as.numeric(values)
    value_range <- range(values, na.rm = TRUE)
    if (!all(is.finite(value_range))) return(rep(NA_real_, length(values)))
    if (diff(value_range) == 0) return(rep((lower + upper) / 2, length(values)))
    lower + (values - value_range[1]) / diff(value_range) * (upper - lower)
  }
  
  # Build a simple spatial grid.
  # Cell size equals neighbourhood_radius.
  # Therefore, any point within the spherical radius of a centre point
  # can only lie in the same or directly adjacent cell along each axis.
  cell_size <- neighbourhood_radius
  
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
  
  # Helper: get candidate indices from nearby grid cells, then apply
  # the spherical Euclidean distance criterion.
  get_sphere_neighbor_indices <- function(i){
    
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
    
    # Exact spherical Euclidean filtering. The surrounding grid cells are
    # only a fast preselection; neighbourhood membership is rotation invariant.
    dx <- x_vec[candidate_idx] - x_vec[i]
    dy <- y_vec[candidate_idx] - y_vec[i]
    dz <- z_vec[candidate_idx] - z_vec[i]
    candidate_idx <- candidate_idx[
      dx^2 + dy^2 + dz^2 <= neighbourhood_radius^2
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
  
  doParallel::registerDoParallel(cores)
  
  # Instead of one foreach task per point, use one foreach task per block.
  # Inside each block, we loop over the center points and accumulate directly.
  k <- NULL  # foreach iteration variable; explicit binding for R CMD check
  normalized_accumulated <- foreach::foreach(k = seq_along(starts),
                                    .combine = combine_norm_results,
                                    .init = list(sum_norm = numeric(n_total),
                                                 count_norm = numeric(n_total))) %dopar% {
                                      
                                      s <- starts[k]
                                      e <- (s + max_rows - 1)
                                      if(e > n_total) e <- n_total
                                      
                                      # Local worker-side accumulators.
                                      # This avoids returning millions of rows.
                                      curr_sum_norm <- numeric(n_total)
                                      curr_count_norm <- numeric(n_total)
                                      
                                      for(i in s:e){
                                        
                                        curr_neighbor_idx <- get_sphere_neighbor_indices(i)
                                        
                                        # Normalize finite values only. Points that never receive a
                                        # valid contribution remain NA in the final result rather than
                                        # being assigned an invented fallback value.
                                        if(length(curr_neighbor_idx) > 1){
                                          curr_values <- col_to_norm_vec[curr_neighbor_idx]
                                          valid_local <- is.finite(curr_values)
                                          if(sum(valid_local) > 1){
                                            valid_idx <- curr_neighbor_idx[valid_local]
                                            valid_values <- curr_values[valid_local]

                                            curr_Q <- stats::quantile(
                                              valid_values,
                                              probs = c(lower_quantile, upper_quantile),
                                              na.rm = TRUE,
                                              names = FALSE
                                            )
                                            curr_values_quantiles <- pmin(
                                              pmax(valid_values, curr_Q[1]),
                                              curr_Q[2]
                                            )
                                            curr_values_normalized <- rescale_local(
                                              curr_values_quantiles, 0, 1
                                            )

                                            curr_sum_norm[valid_idx] <-
                                              curr_sum_norm[valid_idx] + curr_values_normalized
                                            curr_count_norm[valid_idx] <-
                                              curr_count_norm[valid_idx] + 1
                                          }
                                        }
                                      }
                                      
                                      list(sum_norm = curr_sum_norm,
                                           count_norm = curr_count_norm)
                                    }
  
  # terminate cluster
  doParallel::stopImplicitCluster()
  
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
  normalized_values <- rep(NA_real_, n_total)
  has_contribution <- normalized_accumulated$count_norm > 0
  normalized_values[has_contribution] <-
    normalized_accumulated$sum_norm[has_contribution] /
    normalized_accumulated$count_norm[has_contribution]

  df_normalized_summarized <- tibble::tibble(
    n = seq_len(n_total),
    local_height_norm = normalized_values
  )
  
  df_fin <- df %>%
    dplyr::left_join(df_normalized_summarized, by = "n") %>%
    dplyr::select(-n)
  
  # -------------------------------------------------------------------------
  # END OPTIMIZED NORMALIZATION CORE
  # -------------------------------------------------------------------------
  
  # Points without a valid normalization contribution remain NA. This exposes
  # sparse/problematic regions instead of silently assigning a median value.
  if (verbose == TRUE) {
    missing_norm <- sum(!is.finite(df_fin$local_height_norm))
    if (missing_norm > 0) {
      cat(missing_norm, "point(s) could not be normalized and remain NA.\n")
    }
  }
  
  if(verbose == TRUE){
    cat("Adding normalized and contrast-enhanced scales...\n")
  }

  local_height_cols_raw_norm <- get_height_colors(heights = df_fin$local_height_norm)

  # Preserve the historical 10^height emphasis, but rescale its theoretical
  # 1--10 range back to 0--1 so thresholds remain intuitive and comparable.
  local_height_norm_contrast <- (10^df_fin$local_height_norm - 1) / 9
  local_height_norm_contrast <- pmin(pmax(local_height_norm_contrast, 0), 1)
  local_height_cols_contrast_norm <- get_height_colors(heights = local_height_norm_contrast)

  df_fin$local_height_norm_col <- local_height_cols_raw_norm
  df_fin$local_height_norm_contrast <- local_height_norm_contrast
  df_fin$local_height_norm_contrast_col <- local_height_cols_contrast_norm
  
  
  
  if(!is.null(plot_file)){
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("Package 'rgl' is required when plot_file is supplied.", call. = FALSE)
    }
    if(verbose == TRUE){
      cat("2D-plotting to", plot_file, "...\n")
    }
    # Compare the input scale with the locally normalized scale and its
    # peak-enhanced 0--1 contrast representation.
    grDevices::pdf(file = plot_file,
        width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    graphics::par(mfrow=c(1,3))
    input_col <- if (column_to_normalize == "local_height" &&
                     "local_height_col" %in% names(df_fin)) {
      df_fin$local_height_col
    } else {
      get_height_colors(df_fin[[column_to_normalize]])
    }
    graphics::plot(df_fin[[column_to_normalize]], col = input_col, pch=16, cex=.2,
         main = "input",
         xlab = "idx",
         ylab = column_to_normalize)

    graphics::plot(df_fin$local_height_norm, col = df_fin$local_height_norm_col, pch=16, cex=.2,
         main = "normalized",
         xlab = "idx",
         ylab = "local_height_norm (0-1)")

    graphics::plot(df_fin$local_height_norm_contrast, col = df_fin$local_height_norm_contrast_col, pch=16, cex=.2,
         main = "normalized contrast",
         xlab = "idx",
         ylab = "local_height_norm_contrast (0-1)")
    graphics::par(mfrow=c(1,1))
    grDevices::dev.off()
    # dev.print(pdf, file = file.path(df_folder, gsub("csv$", "pdf", curr_filename_out)),
    #           width = (21.0-4)/2.54, height = (29.7-4)/2.54/4)
    
    plot_file_3D <- gsub("normalized\\.pdf", "normalized.png", plot_file)
    
    if(verbose == TRUE){
      print(paste0("Saving 3D comparison of original and normalized files to ", plot_file_3D, "..."))
    }
    
    # plot eye in 'SEM colors'
    rgl::close3d()
    rgl::plot3d(df_fin %>%
             dplyr::select(x,y,z),
           col = df_fin %>%
             dplyr::pull("local_height_col"),
           aspect = "iso",
           size=3)
    
    rgl::points3d(df_fin %>%
               dplyr::select(x,y,z) %>%
               dplyr::mutate(x = x+max(x)+0.5*diff(range(x))),
             col = df_fin %>%
               dplyr::pull("local_height_contrast_col"),
             aspect = "iso",
             size=3)
    
    
    rgl::points3d(df_fin %>%
               dplyr::select(x,y,z) %>%
               dplyr::mutate(z = z-max(z)-0.5*diff(range(z))),
             col = df_fin %>%
               dplyr::pull("local_height_norm_col"),
             aspect = "iso",
             size=3)
    
    rgl::points3d(df_fin %>%
               dplyr::select(x,y,z) %>%
               dplyr::mutate(x = x+max(x)+0.5*diff(range(x)),
                      z = z-max(z)-0.5*diff(range(z))),
             col = df_fin %>% dplyr::pull("local_height_norm_contrast_col"),
             aspect = "iso",
             size=3)
    
    # # View along the X-axis
    # rgl::view3d(userMatrix = rgl::rotationMatrix(pi/2, 0, -1, 0))
    
    rgl::par3d(windowRect = c(20, 30, 800, 800))
    # rgl::par3d(userMatrix = rotate3d(rgl::par3d("userMatrix"), angle2, 0, 0, -1))
    
    # remove bounding box
    # bbox3d(alpha = 0.0, xlab="NULL")
    
    # rotate view roughly to look at eye
    # Step 1: Calculate the mean direction
    mean_vector <- colMeans(df_fin %>%
                              dplyr::select(norm.x, norm.y, norm.z))
    mean_direction <- mean_vector / sqrt(sum(mean_vector^2))  # Normalize
    
    # Step 2: Calculate the opposite direction
    opposite_direction <- -mean_direction
    
    # # Highlight the mean vector and its opposite
    # arrow3d(rep(0, 3), mean_direction*200, col = "red", lwd = 3)
    # arrow3d(rep(0, 3), opposite_direction*200, col = "green", lwd = 3)
    
    # Step 4: Rotate the view to look in the opposite direction of the mean vector
    rgl::view3d(userMatrix = rgl::rotationMatrix(acos(-1),  # 180 degrees rotation
                                       opposite_direction[1],
                                       opposite_direction[2],
                                       opposite_direction[3]))
    
    # # rotate view
    # angle1 <- 45 * (pi/180)   # 45 degrees in radians
    # angle2 <- 60 * (pi/180)   # 60 degrees in radians
    # rgl::par3d(userMatrix = rotate3d(rgl::par3d("userMatrix"), angle1, 0, 1, 0))
    
    rgl::rgl.snapshot(plot_file_3D)
    
    if(plot_results != TRUE){
      rgl::close3d()
    }
  }
  
  if(verbose == TRUE){
    cat("Normalization done!\n")
  }
  
  return(df_fin)
}