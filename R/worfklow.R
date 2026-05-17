#' Convert STL to tibble
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
convert_STL_to_tibble <- function(file_name,
                                  stl_folder,
                                  plot_results = TRUE,
                                  verbose = TRUE){
  # # testing
  # file_name = file_name
  # stl_folder = stl_folder
  # triangle_centers_and_normals_folder = triangle_centers_and_normals_folder
  # plot_results = TRUE
  # write_results_to_file = TRUE
  # verbose = TRUE
  
  if(!grepl("_surface\\.stl", basename(file_name))) stop("File name (", file_name, ") must end with \"_surface.stl\"")
  # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  # fast csv wrighting and reading
  require(readr)
  # dplyr for tibble handling and pipe
  require(dplyr)
  
  # convert mesh data
  
  # define file to work with
  base_file <- gsub("_surface.stl", "", basename(file_name))
  
  if(verbose == TRUE) cat(base_file)
  
  
  # Import and Clean STL and export to csv with coordinates and normals --------
  df <- STL_triangles(file_name = file_name,
                      verbose = verbose)
  
  # multiply by 1000 to account for scaling in Blender
  df <- df %>%
    mutate(x = x * 1000,
           y = y * 1000,
           z = z * 1000,
           norm.x = norm.x * 1000,
           norm.y = norm.y * 1000,
           norm.z = norm.z * 1000)
  
  if(plot_results == TRUE){
    # plot imported eye data
    plot3d(df %>% 
             select(x,y,z),
           size = 5, 
           aspect = "iso")
    
    # draw vectors on eye to see if they point in the right directions
    vec.mult <- .01
    
    x_start <- df %>% 
      select(x)
    y_start <- df %>% 
      select(y)
    z_start <- df %>% 
      select(z)
    
    x_end <- df %>% 
      pull(x) + df %>% 
      pull(norm.x) * vec.mult
    y_end <- df %>% 
      pull(y) + df %>% 
      pull(norm.y) * vec.mult
    z_end <- df %>% 
      pull(z) + df %>% 
      pull(norm.z) * vec.mult
    
    
    # Combine start and end points into a single vector
    segments <- cbind(x_start, y_start, z_start, x_end, y_end, z_end)
    
    # plot corneal projection lines
    segments3d(x=as.vector(t(segments[,c(1,4)])),
               y=as.vector(t(segments[,c(2,5)])),
               z=as.vector(t(segments[,c(3,6)])),
               col = "purple")
  }
  
  if(verbose == TRUE) cat("All done!")
  return(df)
}


#' Find local heights
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'

find_local_heights <- function(df,
                               facet_estimate = 14,
                               cores,
                               plot_results = FALSE,
                               plot_file = NULL,
                               verbose = TRUE,
                               invert = FALSE){
  # # testing
  # facet_estimate = 14
  # cores = 12
  # plot_results = TRUE
  # plot_file = gsub("csv$", "pdf", file_name_out)
  # verbose = TRUE
  # invert = TRUE
  
  # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  # dplyr for tibble handling and pipe
  require(dplyr)
  
  if (verbose == TRUE) cat("Adding local heights according to facet size estimate (", facet_estimate, ")...\n")
  
  # # load csv file
  # df <- read_csv(file_name,
  #                                 show_col_types = FALSE,
  #                                 progress  = FALSE)
  # tri_centers_normals <- df
  
  # add local height to df
  #   This is a multi-threaded but may still take a while. Define number of cores to suit your system (cores = n).
  
  # if(is.null(plot_file)){
  #   df <- calculate_local_heights(df = df,
  #                                 search_diam = facet_estimate*3,
  #                                 cores = cores,
  #                                 verbose = verbose)
  # } else{
    df <- calculate_local_heights(df = df,
                                  search_diam = facet_estimate*3,
                                  cores = cores,
                                  plot_file = plot_file,
                                  verbose = verbose,
                                  invert = invert)
  # }
  
  
  
  # plot eye in 'SEM colors'
  if(plot_results == TRUE){
    plot3d(df %>% 
             select(x,y,z), 
           col = df %>% 
             pull(local_height_log_col), 
           aspect = "iso",
           size=4)
  }
  
  # # write df as csv
  # if(write_results_to_file == TRUE){
  #   write_csv(df,
  #             file_name_out,
  #             progress = FALSE)
  # }
  # if (verbose == TRUE) cat("Local Heights calculated for CV", curr_CV, ".\n")
  # if (verbose == TRUE) cat("*********************************\n")
  
  if (verbose == TRUE) cat("Normalizing local heights.\n")
  df_norm <- normalize_local_heights(df = df,
                                     normalize_diam = facet_estimate,
                                     column_to_normalize = "local_height_log", # "local_height" "local_height_log"
                                     cores = cores,
                                     # plot_file = file.path(df_normalized_folder,
                                     #                       gsub("csv$", "pdf", basename(curr_filename_out))),
                                     verbose = verbose)
  
  if (verbose == TRUE) cat("All done!\n")
  return(df_norm)
}


#' Filter high cuticle vertices by thresholding
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
threshold_high_points <- function(df,
                                  thresholds_file,
                                  column1 = "x",
                                  column2 = "y",
                                  min_threshold = 2,
                                  max_treshold = 8,
                                  final_threshold = 5,
                                  trials = 9,
                                  cores = 18,
                                  plot_results = TRUE,
                                  plot_file,
                                  verbose = FALSE){
  
  # # testing
  # df = local_heights
  # thresholds_file <- "./data/thresholds.log"
  # column1 = "x"
  # column2 = "y"
  # min_threshold = 2
  # max_treshold = 8
  # final_threshold = 5
  # trials = 9
  # cores = 18
  # plot_results = TRUE
  # plot_file = gsub("csv$", "pdf", file_name_out)
  # verbose = TRUE
  
  # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  # dplyr for tibble handling and pipe
  require(dplyr)
  
  # check if a csv file exists to sore all threshold values in a log file
  
  if(!file.exists(thresholds_file)){
    if(verbose == TRUE) cat("Creating ", thresholds_file, "...\n")
    thresholds <- tibble(CV = character(),
                         eye = numeric(),
                         scale = character(),
                         threshold = numeric())
    
    # create empty threshold file
    write_csv(thresholds,
              file = thresholds_file,
              progress = FALSE)
  } else{
    if(verbose == TRUE) cat("Loading ", thresholds_file, "...\n")
    thresholds <- read_csv(file = thresholds_file,
                           show_col_types = FALSE,
                           col_types = "cdcd",
                           progress  = FALSE)
  }
  
  # # load csv file
  # local_heights <- read_csv(file.path(file_name),
  #                           show_col_types = FALSE,
  #                           progress  = FALSE)
  
  # # define new file name independent of raw or normalized input file
  # file_name_out <- gsub("_df.csv", "_df.csv", basename(file_name))
  # file_name_out <- gsub("_df_normalized.csv", "_df.csv", file_name_out)
  
  # # choose from raw values ("local_height") 
  # #   or log-transformed, quantile-filtered values ("local_height_log"), local_height_log_norm
  
  
  # # plot eye in 'SEM colors'
  # if(plot_results == TRUE){
  #   # get columns with highest ranges translations
  #   column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
  #   names(column_ranges) <- c("x","y","z")
  #   column_trans <- names(sort(column_ranges))[2]
  #   plot3d(df %>% 
  #            select(x,y,z), 
  #          col = df %>% 
  #            pull(local_height_col), # gsub("_height", "_height_col", height_column)
  #          aspect = "iso",
  #          size=5)
  #   points3d(df %>% 
  #              select(x,y,z) %>% 
  #              mutate(!!as.symbol(column_trans) := 
  #                       !!as.symbol(column_trans) +
  #                       max(!!as.symbol(column_trans)) +  
  #                       0.2 * diff(range(!!as.symbol(column_trans)))), 
  #            col = df %>% 
  #              pull(local_height_log_col), # gsub("_height", "_height_col", height_column)
  #            aspect = "iso",
  #            size=5)
  #   points3d(df %>% 
  #              select(x,y,z) %>% 
  #              mutate(!!as.symbol(column_trans) := 
  #                       !!as.symbol(column_trans) +
  #                       2 * max(!!as.symbol(column_trans)) +  
  #                       2 * 0.2 * diff(range(!!as.symbol(column_trans)))), 
  #            col = df %>% 
  #              pull(local_height_log_norm_col), # gsub("_height", "_height_col", height_column)
  #            aspect = "iso",
  #            size=5)
  # }
  
  height_column <- "local_height_log_norm" # local_height_log local_height_log_norm
  
  # # show range of height values as orientation for threshoöd selection
  # range(df %>% 
  #         pull(!!as.symbol(height_column)))
  
  # print(file_name)
  if(verbose == TRUE) cat(height_column, "\n")
  
  
  # # find threshold manually
  find_threshold(df = df, # %>% filter(row_number() %% 2 == 0),
                 height_column = height_column,
                 column1 = column1,
                 column2 = column2,
                 min_threshold = min_threshold,
                 max_treshold = max_treshold,
                 trials = trials,
                 plot_file = plot_file)
  
  
  
  
  if(plot_results == TRUE){
    # stop here for manual input of curr_threshold
    # re-run the code from here to # ******************** with changing thresholds
    # until you are satisfied,
    # then continue with the code lines below to save the threshold value
    
    # choose good threshold from plot
    curr_threshold = final_threshold
    
    # filter data according to threshold
    df_thresh <- df %>% 
      filter(!!as.symbol(height_column) >= curr_threshold) %>% 
      select(x,y,z)
    
    if(nrow(df_thresh) > 65536){
      warning("Number (", nrow(df_thresh), ") should not exceed 65536.\nRemove ", nrow(df)-65536, " points or split data.")
    } else{
      if(verbose == TRUE) cat(nrow(df_thresh), " points.\n")
    }
    
    # plot eye in 'SEM colors'
    plot3d(df %>% 
             select(x,y,z), 
           col = df %>% 
             pull(local_height_log_norm_col), 
           aspect = "iso",
           size=10)
    
    # add clusters to plot
    spheres3d(df_thresh, #%>% 
              # slice(110000:111000), 
              col = "orange", 
              radius = 2.5)
  }
  # ********************
  
  # add curr_threshold to thresholds
  thresholds <- thresholds %>% 
    add_row(CV = curr_CV,
            eye = as.numeric(curr_eye),
            scale = height_column,
            threshold = curr_threshold)
  
  # write df with local heights into
  write_csv(thresholds, 
            thresholds_file,
            progress = FALSE)
  
  if(verbose == TRUE) cat("Thresholding done!\n")
  return(df_thresh)
}


#' Find facet poisition candidates
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
#' Find facet positions by thresholding and local condensation
#'
#' This workflow helper first thresholds the local-height table and then runs
#' local mode condensation on the surviving points to obtain one facet
#' candidate per condensed point group.
#'
#' @param df Optional pre-thresholded input table. If supplied and it already
#'   contains `x`, `y`, `z`, and `height_value`, it is used directly.
#' @param local_heights Local-height table used when `df` does not already
#'   contain thresholded points with a `height_value` column.
#' @param local_height_threshold Numeric threshold applied to `height_column`.
#' @param height_column Character string naming the local-height column in
#'   `local_heights`.
#' @param neighbour_radius Local neighbourhood radius used during condensation.
#' @param merge_radius Radius used to merge converged points into final
#'   candidate groups. Defaults to `neighbour_radius / 2`.
#' @param weight_exponent Exponent applied to `height_value` during weighting.
#' @param max_iterations Maximum number of condensation iterations.
#' @param step_size Fraction of each shift toward the local weighted centre.
#' @param tolerance Early-stop tolerance for the maximum per-iteration shift.
#' @param min_cluster_size Minimum cluster size retained as a candidate.
#' @param select_point How to choose the final candidate point within a
#'   converged group. Passed to `find_facet_candidates_condensed()`.
#' @param return_details Logical. If `TRUE`, return the full details list from
#'   `find_facet_candidates_condensed()`. Otherwise return only the candidates
#'   table.
#' @param verbose Logical. If `TRUE`, print short progress information.
#' @param ... Currently unused. Kept so old project scripts can be adapted more
#'   gently without failing on extra arguments immediately.
#'
#' @return Data frame of facet candidates, or the full details list if
#'   `return_details = TRUE`.
#'
#' @export
find_facet_positions <- function(df = NULL,
                                 local_heights = NULL,
                                 local_height_threshold = 2.5,
                                 height_column = "local_height_log",
                                 neighbour_radius,
                                 merge_radius = neighbour_radius / 2,
                                 weight_exponent = 1,
                                 max_iterations = 10,
                                 step_size = 1,
                                 tolerance = neighbour_radius * 1e-3,
                                 min_cluster_size = 1,
                                 select_point = c("nearest_mode", "max_height"),
                                 return_details = FALSE,
                                 verbose = FALSE,
                                 ...){

  # # testing
  # local_height_threshold = 2.5
  # height_column = "local_height_log"
  # neighbour_radius = 15
  # merge_radius = 7.5
  # verbose = TRUE

  select_point <- match.arg(select_point)

  has_thresholded_input <- is.data.frame(df) &&
    all(c("x", "y", "z", "height_value") %in% colnames(df))

  if(has_thresholded_input){
    thresholded_points <- df

    # Keep source_index available for traceability if not already present.
    if(!"source_index" %in% colnames(thresholded_points)){
      thresholded_points$source_index <- seq_len(nrow(thresholded_points))
    }

    thresholded_points <- thresholded_points %>%
      dplyr::select(source_index, x, y, z, height_value)

  } else {

    if(is.null(local_heights)){
      stop(
        "Either provide df with x/y/z/height_value or provide local_heights ",
        "together with height_column.",
        call. = FALSE
      )
    }

    thresholded_points <- threshold_local_heights(
      df = local_heights,
      local_height_threshold = local_height_threshold,
      height_column = height_column,
      verbose = verbose
    )
  }

  facet_candidates <- find_facet_candidates_condensed(
    df = thresholded_points,
    coord_cols = c("x", "y", "z"),
    height_col = "height_value",
    neighbour_radius = neighbour_radius,
    merge_radius = merge_radius,
    weight_exponent = weight_exponent,
    max_iterations = max_iterations,
    step_size = step_size,
    tolerance = tolerance,
    min_cluster_size = min_cluster_size,
    select_point = select_point,
    return_details = return_details,
    verbose = verbose
  )

  return(facet_candidates)
}

#' Calculate optic parameters for eye
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'

get_optic_parameters <- function(df,
                                 plot_file = NULL,
                                 edge_tol = 0.5,
                                 cores = 1,
                                 plot_results = FALSE,
                                 verbose = FALSE,
                                 version = "fast",
                                 facet_size = 14){
  
  # # testing
  # df = LMs_facets_combined
  # plot_file = file.path(facet_infos_folder,
  #                       gsub("_surface.stl", "_neighbour_and_size_data.pdf", basename(file_name)))
  # edge_tol = 0.8
  # cores = 18
  # plot_results = TRUE
  # verbose = TRUE
  # facet_size = facet_estimate
  # version = "fast"
  # version = "slow"
  
  # get data of facets
  curr_facets <- df %>% 
    filter(type == "facet")
  
  # get data of LMs
  curr_LMs <- df %>% 
    filter(type == "LM")
  
  
  if(version == "fast"){
    df_neighbours <- find_neighbours(df = curr_facets,
                                     edge_tol = edge_tol)
  } else{
    df_neighbours_raw <- find_neighbours_deprecated(df = curr_facets,
                                                    facet_size = facet_estimate,
                                                    neighbour_threshold = 1.5,
                                                    cores = cores,
                                                    plot_results = plot_results,
                                                    plot_file = NULL,
                                                    verbose = verbose)
    
    
    
    # # create viridis colors
    # color_option = "D"
    # pal <- viridisLite::viridis(6, option = color_option)
    # 
    # #     add data - replace size by size avg
    # df_neighbours <- curr_facets %>% 
    #   left_join(df_neighbours_raw %>% 
    #               select(-c(size, size_avg)) ,
    #             by = "ID") %>% 
    #   mutate(number_of_neighs_cols = pal[pmin(pmax(number.of.neighbours, 1), 6)])
    
    df_neighbours <- curr_facets %>%
      left_join(df_neighbours_raw, # %>%
                # select(-c(size, size_avg)) ,
                by = "ID") #%>%
    # mutate(number_of_neighs_cols = pal[pmin(pmax(number.of.neighbours, 1), 6)])
  }
  
  # add results to tibble
  df_w_neighbours <- curr_facets %>% 
    left_join(df_neighbours %>% 
                select(-c(x,y,z,type)), 
              by = "ID")
  
  facet_sizes <- calculate_facet_size(df_w_neighbours)
  
  
  
  df_w_sizes_raw <- left_join(df_w_neighbours,
                              facet_sizes %>% 
                                select(-n_used),
                              by = "ID")
  
  #   go on with avg sizes
  if(verbose == TRUE) cat("Using average size...\n")
  df_w_sizes <- df_w_sizes_raw %>% 
    select(-size) %>% 
    rename(size = size_avg)
  
  # #   check:
  # plot3d(df_w_sizes %>%
  #          select(x,y,z),
  #        col = df_w_sizes$number_of_neighs_cols,
  #        type="s",
  #        radius = df_w_sizes$size,
  #        aspect = "iso")
  # 
  # 
  # texts3d(df_w_sizes  %>%
  #           select(x,y,z),
  #         texts = df_w_sizes  %>%
  #           pull(ID),
  #         pos=1,
  #         cex = .7)
  
  
  # calculate facet normals according to their neighbours
  df_normals <- get_facet_normals(df = df_w_sizes,
                                  cores = cores,
                                  plot_file = gsub("_neighbour_and_size_data", 
                                                   "_normal_data", 
                                                   plot_file),
                                  plot_results = TRUE,
                                  verbose = TRUE)
  
  
  
  
  # add info to df
  df_w_normals <- df_w_sizes %>% 
    left_join(df_normals, 
              by="ID")
  
  
  # #   check:
  # plot3d(df_w_normals %>%
  #          select(x,y,z),
  #        col = df_w_normals$number_of_neighs_cols,
  #        type="s",
  #        radius = df_w_normals$size,
  #        aspect = "iso")
  # 
  # 
  # make_segments(df = df_w_normals,
  #               start_colums = c(2,3,4),
  #               end_colums = c(10,11,12),
  #               vector_length_multipler = 500)
  
  # calculateIF angle, P, CPD
  optic_parameters <- get_optic_properties(df = df_w_normals,
                                           cores = cores,
                                           plot_results = TRUE,
                                           plot_file = gsub("_neighbour_and_size_data", 
                                                            "_optics_parameters", 
                                                            plot_file),
                                           verbose = TRUE)
  
  
  
  # add results to tibble
  df_w_optic_parameters <- df_w_normals %>% 
    left_join(optic_parameters, by = "ID")
  
  
  # #   check:
  # plot3d(df_w_optic_parameters %>%
  #          select(x,y,z),
  #        col = df_w_optic_parameters$number_of_neighs_cols,
  #        type="s",
  #        radius = df_w_optic_parameters$size,
  #        aspect = "iso")
  # 
  # 
  # texts3d(df_w_optic_parameters  %>%
  #           select(x,y,z),
  #         texts = df_w_optic_parameters  %>%
  #           pull(ID),
  #         pos=1,
  #         cex = .7)
  
  # define viridis colours for the different parameters
  cols_to_use <- viridis(n=100, begin = 0, end = 1)
  
  df_w_plot_cols <- df_w_optic_parameters %>% 
    arrange(size) %>% 
    mutate(size_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                               arrange(size) %>% 
                                               pull(size),
                                             cols_to_use)) %>%
    arrange(P) %>% 
    mutate(P_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                            arrange(P) %>% 
                                            pull(P),
                                          cols_to_use)) %>%
    arrange(v) %>% 
    mutate(v_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                            arrange(v) %>% 
                                            pull(v),
                                          cols_to_use)) %>%
    arrange(CPD) %>% 
    mutate(CPD_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                              arrange(CPD) %>% 
                                              pull(CPD),
                                            cols_to_use)) %>%
    arrange(delta_phi.deg) %>% 
    mutate(delta_phi.deg_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                                        arrange(delta_phi.deg) %>% 
                                                        pull(delta_phi.deg),
                                                      cols_to_use)) %>%
    arrange(number.of.neighbours) %>% 
    mutate(number_of_neighs_cols = continuous_color_ramp(df_w_optic_parameters %>% 
                                                           arrange(number.of.neighbours) %>% 
                                                           pull(number.of.neighbours),
                                                         cols_to_use))
  
  
  # #   check:
  # plot3d(df_w_plot_cols %>%
  #          select(x,y,z),
  #        col = df_w_plot_cols$size_cols,
  #        # col = df_w_plot_cols$number_of_neighs_cols,
  #        # col = df_w_plot_cols$delta_phi.deg_cols,
  #        # col = df_w_plot_cols$CPD_cols,
  #        type="s",
  #        radius = df_w_plot_cols$size,
  #        aspect = "iso")
  # 
  # 
  # texts3d(df_w_plot_cols  %>%
  #           select(x,y,z),
  #         texts = df_w_plot_cols  %>%
  #           pull(ID),
  #         pos=1,
  #         cex = .7)
  
  return(df_w_plot_cols)
}


#' Normalize data of all eyes
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
normalize_eye_data <- function(df,
                               info_table = CompVisTab, 
                               quantile_value = 0.01,
                               plot_results = FALSE){
  
  # # testing
  # df = df_all
  # info_table = CompVisTab
  # quantile_value = 0.01
  # plot_results = TRUE
  
  # for viridis colours
  require(viridis)
  
  # fast csv writing and reading
  require(readr)
  
  # load tidyverse for its various conveniences
  require(dplyr)
  rename <- dplyr::rename
  
  # add species to df_all
  df_all <- df_all %>% 
    left_join(info_table %>% 
                mutate(taxon = paste(CV, genus, species)) %>% # CV, subfamily, tribe, 
                select(CV, taxon) %>% # , subfamily, tribe
                distinct(), by = "CV")
  
  # define analysis variables
  directions <- c("latitude", "longitude")
  analysis_variables <- c("size", "delta_phi.deg", "P", "CPD")
  analysis_variables_indices <- paste0(0, 1:4)
  
  
  # colorize all corrected data values with table-wide LUT
  cols_to_use <- viridis(n=200, begin = 0, end = 1)
  
  # first: remove outliers
  df_all_corr <- NULL
  
  curr_CV <- unique(df_all$CV)[1]
  for(curr_CV in unique(df_all$CV)){
    curr_df_all <- df_all %>% 
      filter(CV == curr_CV)
    
    curr_var <- analysis_variables[1]
    for(curr_var in analysis_variables){
      curr_values <- curr_df_all %>%  #unlist(df_all[,curr_var_col]) 
        pull(get(!!!curr_var))
      curr_NA_vals <- which(is.na(curr_values))
      curr_values[curr_NA_vals] <- mean(curr_values, na.rm = TRUE)
      
      # find outliers and replace them by quantile values
      curr_quantiles <- quantile(curr_values, c(quantile_value, 1-quantile_value))
      curr_values[curr_values<curr_quantiles[1]] <- curr_quantiles[1]
      curr_values[curr_values>curr_quantiles[2]] <- curr_quantiles[2]
      
      curr_df_all <- curr_df_all %>% 
        mutate(new_col = curr_values)
      colnames(curr_df_all)[ncol(curr_df_all)] <- paste0(curr_var, "_corr")
      
      curr_df_all <- curr_df_all %>% 
        mutate(new_col = continuous_color_ramp(curr_values,
                                               cols_to_use))
      colnames(curr_df_all)[ncol(curr_df_all)] <- paste0(curr_var, "_corr_cols")
      
      curr_df_all[curr_NA_vals,ncol(curr_df_all)] <- NA
    }
    df_all_corr <- rbind(df_all_corr, curr_df_all)
  }
  
  # re-define variables to use outlier-free data
  analysis_variables <- c("size_corr", "delta_phi.deg_corr", "P_corr", "CPD_corr")
  curr_var <- analysis_variables[1]
  for(curr_var in analysis_variables){
    curr_values <- df_all_corr %>%
      pull(get(!!!curr_var))
    curr_NA_vals <- which(is.na(curr_values))
    curr_values[curr_NA_vals] <- mean(curr_values, na.rm = TRUE)
    
    df_all_corr <- df_all_corr %>% 
      mutate(new_col = continuous_color_ramp(curr_values,
                                             cols_to_use))
    colnames(df_all_corr)[ncol(df_all_corr)] <- paste0(curr_var, "_cols_all")
    
    curr_values[curr_NA_vals] <- NA
  }
  
  
  # plot to compare raw, outlier-free, and normalized data
  if(plot_results == TRUE){
    par(mfrow=c(3,1))
    plot(df_all_corr %>%
           filter(!is.na(size)) %>%
           pull(x), col = df_all_corr$CPD_cols, pch=16)
    plot(df_all_corr %>%
           filter(!is.na(size)) %>%
           pull(x), col = df_all_corr$CPD_corr_cols, pch=16)
    plot(df_all_corr %>%
           filter(!is.na(size)) %>%
           pull(x), col = df_all_corr$CPD_corr_cols_all, pch=16)
    par(mfrow=c(1,1))
  }
  
  # add variables to info_table
  # mean facet size
  info_table <- info_table %>% 
    left_join(df_all_corr %>% 
                group_by(CV) %>% 
                summarise(mean_facet_size = mean(size_corr, na.rm = TRUE)),
              by = "CV")
  
  # mean IF angle
  info_table <- info_table %>% 
    left_join(df_all_corr %>% 
                group_by(CV) %>% 
                summarise(mean_IF_angle = mean(delta_phi.deg_corr, na.rm = TRUE)),
              by = "CV")
  # mean P
  info_table <- info_table %>% 
    left_join(df_all_corr %>% 
                group_by(CV) %>% 
                summarise(mean_P = median(P_corr, na.rm = TRUE)), # here!: mean-median
              by = "CV")
  
  # mean V
  info_table <- info_table %>% 
    left_join(df_all_corr %>% 
                group_by(CV) %>% 
                summarise(mean_CPD = median(CPD_corr, na.rm = TRUE)), # here!: mean-median
              by = "CV")
  
  # total surface area of eye using facet sizes
  info_table <- info_table %>% 
    left_join(df_all_corr %>% 
                group_by(CV) %>% 
                mutate(eye_surface_area = sum(calculate_hexagon_area(size_corr))) %>% 
                select(CV, eye_surface_area) %>% 
                distinct(),
              by = "CV") 
  
  
  # add eye surface area to df_all_corr
  df_all_corr <- df_all_corr %>% 
    left_join(info_table %>% 
                distinct(CV, .keep_all = TRUE) %>%
                select(CV, eye_surface_area),
              by = "CV")
  
  # create new df from existing one. xxx: rename info_table even earlier, after initial loading
  info_table_filtered <- info_table
  
  
  # data cleaning -----------------------------------------------------------
  df_all_corr_export <- df_all_corr %>% 
    select(-c(filename ,
              size, delta_phi.rad, delta_phi.deg, P, v, CPD, 
              size_cols, P_cols, v_cols, CPD_cols, delta_phi.deg_cols, delta_phi.deg))
  
  return(list(info_table_filtered,
              df_all_corr_export))
}




#' Calculate Corneal Projections
#'
#' XYZ: create description and param defs
#'
#' @param df A tibble containing coordinates and IDs.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
calculate_corneal_projections <- function(df,
                                          cp_diam_cm,
                                          center = NULL,
                                          verbose = FALSE){
  # #   testing
  # df = df_all_normalized
  # cp_diam_cm = 50
  # center = head_center
  # verbose = TRUE
  
  # General variable definitions
  cp_area_um2 = 4*pi*((cp_diam_cm/2)*1000)^2
  
  
  # corneal projections, FOV and Voronoi -----------------------------------------------------------------
  
  # add latitude and longitude columns
  if(all(grepl("latitude", colnames(df)) == FALSE) == TRUE){
    df$latitude <- NA
    df$longitude <- NA
  }
  
  # add corn. projection columns
  if(all(grepl("corn.proj.x", colnames(df)) == FALSE) == TRUE){
    df$corn.proj.x <- NA
    df$corn.proj.y <- NA
    df$corn.proj.z <- NA
  }
  
  # df %>% 
  #   distinct(CV) %>% 
  #   print(n=999)
  
  k=1
  for(k in 1:length(unique(df$CV))){ # 1:length(unique(df$CV))
    # get current CV
    curr_CV <- unique(df$CV)[k]
    if(verbose == TRUE) cat(paste0("CV",curr_CV), "\n")
    
    # get data for each facet of current CV
    curr_df <- df %>% 
      filter(CV == curr_CV,
             type == "facet") %>% 
      select(-c(latitude, longitude))# %>%
    
    # calculate center of left eye
    center_point_L <- c(mean(curr_df$x),
                        mean(curr_df$y),
                        mean(curr_df$z))
    
    
    # define center of head - this is always the center od the global coord system
    center_point_head <- c(0, 0, 0)
    
    # center between eyes; eye distance radius value
    # define center of corneal projection sphere
    #   this is exaclty between both eyes, so x=0:
    
    if(is.null(center)){
      sphere.c <- center_point_L
      sphere.c[1] <- 0
    } else{
      sphere.c <- center
    }
    # calculate sphere radius in um
    sphere.r = cp_diam_cm/2*1000 # abs(center_point_L[1]) + 
    
    # plot the eye with IF angle colors
    plot3d(curr_df %>% 
             select(x,y,z), 
           radius = curr_df$size_corr, 
           col = curr_df$CPD_corr_cols,
           type="s", 
           label = T, 
           add = F, 
           aspect = "iso")
    
    # add the landmarks in blue
    spheres3d(df %>%
                filter(CV == curr_CV,
                       type == "LM") %>%
                select(x,y,z),
              radius = 2*mean(curr_df$size_corr), 
              col="blue")
    
    # add corneal projection sphere
    spheres3d(sphere.c, 
              radius = sphere.r, 
              col=rgb(0,0,0.3), alpha=0.05)
    
    
    # calculate normal intersections with sphere left
    if(verbose == TRUE) cat("Calculating corneal projection intersections...\n")
    
    curr_df$corn.proj.x <- NA
    curr_df$corn.proj.y <- NA
    curr_df$corn.proj.z <- NA
    for(i in 1:nrow(curr_df)){
      point = c(curr_df$x[i], curr_df$y[i], curr_df$z[i])
      vector = c(curr_df$norm.x[i], curr_df$norm.y[i], curr_df$norm.z[i])
      intersection <- vector.sphere.intersect(point, vector, sphere.c, sphere.r)[2, ]
      curr_df$corn.proj.x[i] <- intersection[1] %>% pull()
      curr_df$corn.proj.y[i] <- intersection[2] %>% pull()
      curr_df$corn.proj.z[i] <- intersection[3] %>% pull()
    }
    
    # add corn. projection coordinates to df
    df$corn.proj.x[df$CV == curr_CV &
                     df$type == "facet"] <- curr_df$corn.proj.x
    df$corn.proj.y[df$CV == curr_CV &
                     df$type == "facet"] <- curr_df$corn.proj.y
    df$corn.proj.z[df$CV == curr_CV &
                     df$type == "facet"] <- curr_df$corn.proj.z
    
    # # draw corn proj vectors
    # vec.mult <- sphere.r # /1000
    # for(l in seq(1, nrow(curr_df), length.out = 100)){
    #   lines3d(x = c(curr_df$x[l], curr_df$x[l] +
    #                   curr_df$norm.x[l]*vec.mult),
    #           y = c(curr_df$y[l], curr_df$y[l] +
    #                   curr_df$norm.y[l]*vec.mult),
    #           z = c(curr_df$z[l], curr_df$z[l] +
    #                   curr_df$norm.z[l]*vec.mult),
    #           col = curr_df$CPD_corr_cols[l])
    # }
    
    # calculate latitudes and longitudes for the curr data
    lat_lon <- convert_to_latlon(x = curr_df %>% 
                                   pull(corn.proj.x),
                                 y = curr_df %>% 
                                   pull(corn.proj.y),
                                 z = curr_df %>% 
                                   pull(corn.proj.z))
    
    # add lat lon data to df
    df$latitude[df$CV == curr_CV &
                  df$type == "facet"] <- lat_lon$latitude
    df$longitude[df$CV == curr_CV &
                   df$type == "facet"] <- lat_lon$longitude
  }
  return(df)
}