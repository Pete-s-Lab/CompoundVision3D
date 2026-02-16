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
                               verbose = TRUE){
  # # testing
  # facet_estimate = 14
  # cores = 12
  # plot_results = TRUE
  # plot_file = gsub("csv$", "pdf", file_name_out)
  # verbose = TRUE
  
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
  
  if(is.null(plot_file)){
    df <- calculate_local_heights(df = df,
                       search_diam = facet_estimate*3,
                       cores = cores,
                       verbose = verbose)
  } else{
    df <- calculate_local_heights(df = df,
                       search_diam = facet_estimate*3,
                       cores = cores,
                       plot_file = plot_file,
                       verbose = verbose)
  }
  
  
  
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
                 trials = 9,
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
find_facet_positions <- function(df,
                                 local_heights,
                                 h_min = 40,
                                 h_max = 16,
                                 h_final = 10,
                                 trials = 9,
                                 column1 = NULL,
                                 column2 = NULL,
                                 cutoffs_file,
                                 plot_results = TRUE,
                                 plot_file,
                                 verbose = FALSE){
  # testing
  # h_min = 20
  # h_max = 8
  # h_final = 16
  # cutoffs_file = "./data/cutoffs.log"
  # plot_results = TRUE
  
  # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  # dplyr for tibble handling and pipe
  require(dplyr)
  
  
  # check if a csv file exists to sore all threshold values in a log file
  if(!file.exists(cutoffs_file)){
    if(verbose == TRUE) cat("Creating ", cutoffs_file, "...\n")
    cutoffs <- tibble(CV = character(),
                      eye = numeric(),
                      cutoff_min = numeric(),
                      cutoff_max = numeric(),
                      cutoff_final = numeric())
    
    # create empty threshold file
    write_csv(cutoffs, file = cutoffs_file,
              progress  = FALSE)
  } else{
    if(verbose == TRUE) cat("Loading ", cutoffs_file, "...\n")
    cutoffs <- read_csv(file = cutoffs_file,
                        show_col_types = FALSE,
                        col_types = "cdddd",
                        progress  = FALSE)
  }
  
  
  # Find facet position candidates
  
  # # load csv file of rough clusters
  # rough_clusters <- read_csv(file_name,
  #                            show_col_types = FALSE)
  
  # # get local heights for plotting
  # curr_df_file <- gsub(df_folder, df_folder, file_name)
  # curr_df_file <- gsub("_df", "_df", curr_df_file)
  # df <- read_csv(curr_df_file,
  #                           show_col_types = FALSE)
  
  
  if(plot_results == TRUE){
    # plot local heights
    plot3d(local_heights %>% 
             select(x,y,z),
           aspect = "iso",
           col = local_heights$local_height_log_col,
           size = 5)
    
    # plot rough clusters
    points3d(df,
             size = 10,
             col = "orange")
  }
  
  # get fine peaks. Make sure your plot device is as large as possible
  facet_positions_auto <- find_facet_canidates(df = df, # %>% 
                                               # slice(1:floor(nrow(.)/2)), 
                                               # slice(ceiling(nrow(.)/2):nrow(.)), 
                                               h_min = h_min,
                                               h_max = h_max,
                                               column1 = NULL,
                                               column2 = NULL,
                                               h_final = h_final,
                                               n_steps = 100,
                                               trials = trials,
                                               plot_file = plot_file,
                                               verbose = verbose)
  
  
  # # in case analysis was split in two halves (done for CV0034)
  # facet_positions_auto_02 <- facet_positions_auto_02 %>% 
  #   mutate(ID = ID+nrow(facet_positions_auto_01))
  # facet_positions_auto <- rbind(facet_positions_auto_01, facet_positions_auto_02)
  
  if(plot_results == TRUE){
    # plot facet position candidates over rough clusters to check again if necessary
    # plot local heights
    plot3d(local_heights %>% 
             select(x,y,z),
           aspect = "iso",
           col = local_heights$local_height_log_col,
           size = 5)
    
    # plot rough clusters
    points3d(df,
             size = 10,
             col = "orange")
    
    # plot facet position candidates
    spheres3d(facet_positions_auto %>% 
                select(x,y,z), 
              col = "red", 
              radius=3)
  }
  
  # add values to cutoffs and save csv
  cutoffs <- cutoffs %>% 
    add_row(CV = curr_CV,
            eye = as.numeric(curr_eye),
            cutoff_min = facet_positions_auto$cutoff_min[1],
            cutoff_max = facet_positions_auto$cutoff_max[1],
            cutoff_final = facet_positions_auto$cutoff_fin[1])
  
  # write cutoffs to csv
  write_csv(cutoffs, 
            cutoffs_file,
            progress = FALSE)
  
  # # write facet_positions_auto
  # write_csv(facet_positions_auto %>% 
  #             select(ID,x,y,z), 
  #           file.path(facet_candidate_folder, file_name_out[1]),
  #           progress = FALSE)
  
  if(verbose == TRUE) cat("All done!\n")
  
  if(verbose == TRUE) cat("Now go to Blender and check facet candidates at", file.path(facet_candidate_folder, file_name_out[1]), "manually.\n")
  
  return(facet_positions_auto %>% 
           select(ID,x,y,z))
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
                                 cores = 18,
                                 plot_results = FALSE,
                                 verbose = FALSE){
  
  # # testing
  # df = LMs_facets_combined
  # plot_file = file.path(facet_infos_folder,
  #                       gsub("_surface.stl", "_neighbour_and_size_data.pdf", basename(file_name)))
  # cores = 18
  # plot_results = TRUE
  # verbose = TRUE
  
  # get data of facets
  curr_facets <- df %>% 
    filter(type == "facet")
  
  # get data of LMs
  curr_LMs <- df %>% 
    filter(type == "LM")
  
  df_neighbours <- find_neighbours(df = curr_facets,
                                   edge_tol = 0.5)
  
  
  # add results to tibble
  df_w_neighbours <- curr_facets %>% 
    left_join(df_neighbours %>% 
                select(-c(x,y,z,type)), 
              by = "ID")
  
  facet_sizes <- calculate_facet_size(df_w_neighbours)
  
  
  
  df_w_sizes <- left_join(df_w_neighbours,
                          facet_sizes %>% 
                            select(-n_used),
                          by = "ID") %>% 
    rename(size = facet_size)
  
  
  # #   check:
  #   plot3d(df_w_sizes %>% 
  #            select(x,y,z),
  #          col = df_w_sizes$number_of_neighs_cols,
  #          type="s",
  #          radius = df_w_sizes$size,
  #          aspect = "iso")
  #   
  #   
  #   texts3d(df_w_sizes  %>% 
  #             select(x,y,z),
  #           texts = df_w_sizes  %>% 
  #             pull(ID),
  #           pos=1,
  #           cex = .7)
  
  
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
  #        # col = df_w_plot_cols$number_of_neighs_cols,
  #        # col = df_w_plot_cols$delta_phi.deg_cols,
  #        col = df_w_plot_cols$CPD_cols,
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
                               plot_results = FALSE){
  
  # # testing
  # df = df_all
  # info_table = info_table
  # plot_results = TRUE
  
  # for viridis colours
  require(viridis)
  
  # fast csv writing and reading
  require(readr)
  
  # load tidyverse for its various conveniences
  require(dplyr)
  rename <- dplyr::rename
  
  # functions
  
  ## -------------------------------------------------------------------------
  ## General variable definitions
  
  ## -------------------------------------------------------------------------
  
  
  
  # add CV number to df_all
  df_all <- df %>% 
    mutate(CV = gsub("^CV(.{4}).*", "\\1", df_all$filename), 
           .before = ID) %>% 
    # arrange(CV, ID)
    arrange(CV, order(gtools::mixedorder(ID)))
  
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
  cols_to_use <- viridis(n=1000, begin = 0, end = 1)
  
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
      curr_quantiles <- quantile(curr_values, c(.05, .95))
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
                                          verbose = FALSE){
  # #   testing
  # df = df_all_normalized
  # cp_diam_cm = 100
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
    sphere.c <- center_point_L
    sphere.c[1] <- 0
    
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