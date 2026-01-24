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
                               curr_facet_estimate = 14,
                               cores,
                               plot_results = FALSE,
                               plot_file = NULL,
                               verbose = TRUE){
  # # testing
  # curr_facet_estimate = 14
  # cores = 12
  # plot_results = TRUE
  # plot_file = gsub("csv$", "pdf", file_name_out)
  # verbose = TRUE
  
  # Dependencies ------------------------------------------------------------
  # 3D plotting
  require(rgl)
  # dplyr for tibble handling and pipe
  require(dplyr)
  
  if (verbose == TRUE) cat("Adding local heights according to facet size estimate (", curr_facet_estimate, ")...\n")
  
  # # load csv file
  # df <- read_csv(file_name,
  #                                 show_col_types = FALSE,
  #                                 progress  = FALSE)
  # tri_centers_normals <- df
  
  # add local height to df
  #   This is a multi-threaded but may still take a while. Define number of cores to suit your system (cores = n).
  
  if(is.null(plot_file)){
  df <- calculate_df(df = df,
                                           search_diam = curr_facet_estimate*3,
                                           cores = cores,
                                           verbose = verbose)
  } else{
    df <- calculate_df(df = df,
                                             search_diam = curr_facet_estimate*3,
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
  df_norm <- normalize_df(df = df,
                                                normalize_diam = curr_facet_estimate,
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
  # thresholds_file <- "./data/thresholds.log"
  # column1 = "x"
  # column2 = "y"
  # min_threshold = 2
  # max_treshold = 8
  # final_threshold = 5
  # cores = 18
  # plot_results = TRUE
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
  
  
  # plot eye in 'SEM colors'
  if(plot_results == TRUE){
    # get columns with highest ranges translations
    column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
    names(column_ranges) <- c("x","y","z")
    column_trans <- names(sort(column_ranges))[2]
    plot3d(df %>% 
             select(x,y,z), 
           col = df %>% 
             pull(local_height_col), # gsub("_height", "_height_col", height_column)
           aspect = "iso",
           size=5)
    points3d(df %>% 
               select(x,y,z) %>% 
               mutate(!!as.symbol(column_trans) := 
                        !!as.symbol(column_trans) +
                        max(!!as.symbol(column_trans)) +  
                        0.2 * diff(range(!!as.symbol(column_trans)))), 
             col = df %>% 
               pull(local_height_log_col), # gsub("_height", "_height_col", height_column)
             aspect = "iso",
             size=5)
    points3d(df %>% 
               select(x,y,z) %>% 
               mutate(!!as.symbol(column_trans) := 
                        !!as.symbol(column_trans) +
                        2 * max(!!as.symbol(column_trans)) +  
                        2 * 0.2 * diff(range(!!as.symbol(column_trans)))), 
             col = df %>% 
               pull(local_height_log_norm_col), # gsub("_height", "_height_col", height_column)
             aspect = "iso",
             size=5)
  }
  
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
    df <- df %>% 
      filter(!!as.symbol(height_column) >= curr_threshold) %>% 
      select(x,y,z)
    
    if(nrow(df) > 65536){
      warning("Number (", nrow(df), ") should not exceed 65536.\nRemove ", nrow(df)-65536, " points or split data.")
    } else{
      if(verbose == TRUE) cat(nrow(df), " points.\n")
    }
    
    # plot eye in 'SEM colors'
    plot3d(df %>% 
             select(x,y,z), 
           col = df %>% 
             pull(local_height_log_col), 
           aspect = "iso",
           size=10)
    
    # add clusters to plot
    spheres3d(df, #%>% 
              # slice(110000:111000), 
              col = "orange", 
              aspect = "iso",
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
  
  if(verbose == TRUE) cat("All done!\n")
  return(df)
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
                                 h_min = 46.127,
                                 h_max = 16.287,
                                 h_final = 23.873,
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

