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





#' Rotate and translate facet positions according to the global coordinate system
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
rotate_and_translate_to_global <- function(df,
                                           facet_estimate = 14,
                                           cores = 18,
                                           plot_results = TRUE,
                                           verbose = TRUE){
  curr_file = file_list_facets_LMs[CV]
  print(curr_file)
  
  curr_CV <- gsub("^CV(\\d+)_.+", "\\1", curr_file)
  curr_eye <- gsub("^.+eye(\\d+)_.+$", "\\1", curr_file)
  
  # get current facet size estimate
  curr_facet_estimate <- eye_infos_df %>%
    filter(CV == curr_CV, eye == curr_eye) %>% 
    pull(facet_size_estimate)
  
  # load all data
  curr_data <- read_csv(file.path(facets_landmarks_combined, curr_file),
                        show_col_types = FALSE)
  
  # get data of facets
  curr_facets <- curr_data %>% 
    filter(type == "facet")
  
  # get data of LMs
  curr_LMs <- curr_data %>% 
    filter(type == "LM")
  
  # get data for left eye
  eye_L <- curr_facets %>% 
    mutate(facet = row_number(),
           ID = paste0("L_", ID)) %>%
    select(c(ncol(.), 1:(ncol(.)-1)))
  
  # # do the next calculations only for left eye and mirror results to right eye
  
  # detach("package:CompoundVision3D", unload=TRUE)
  # install.packages("CompoundVision3D")
  library(CompoundVision3D)
  
  neighbours <- find_neighbours(df = eye_L,
                                facet_size <- curr_facet_estimate,
                                cores = 12,
                                plot_file = file.path(facet_infos, 
                                                      gsub("_facets_LMs_combined.csv", "_neighbour_and_size_data.pdf", curr_file)),
                                verbose = TRUE)
  
  # add results to tibble
  eye_L <- eye_L %>% 
    left_join(neighbours, by = "facet")
  
  # create color range for plotting
  cols.no.neighbours <- continuous_color_ramp(eye_L$number.of.neighbours,
                                              viridis(n=6))
  # eye_L$cols_no_neighbours <- cols.no.neighbours
  
  cols.size <- continuous_color_ramp(eye_L$size, 
                                     viridis((n=9)))
  
  # plot results
  plot3d(eye_L %>%
           select(x,y,z),
         radius = eye_L$size,
         col = cols.no.neighbours,
         type="s",
         label = TRUE,
         aspect = "iso") # , alpha = 0.2
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  plot3d(eye_L %>% 
           select(x,y,z),
         radius =  eye_L$size, # eye_L$size, # 10, #
         type="s",
         label = TRUE,
         aspect = "iso",
         col = cols.size)
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  
  # print("Saving rotated data png...")
  # # rgl.snapshot(file.path("data/1_pre_STLs/CV0002_Adult_Drosohila_1.2um_zsc50/CV0002_Adult_Drosohila_1.2um_zsc50_eye2/stl/imgs", paste0("neighbors_", curr_CV, ".png")))
  
  
  
  # calculate facet normals according to their neighbors
  eye_L_normals <- get_facet_normals(df = eye_L,
                                     cores = 12,
                                     plot_file = file.path(facet_infos,
                                                           gsub("_facets_LMs_combined.csv", "_normal_data.pdf", curr_file)),
                                     verbose = TRUE)
  
  
  
  # add info to eye_L
  eye_L <- eye_L %>% 
    left_join(eye_L_normals, 
              by="facet")
  
  
  # calculateIF angle, P, CPD
  optic_parameters <- get_optic_properties(df = eye_L,
                                           cores = 12,
                                           plot_file = file.path(facet_infos,
                                                                 gsub("_facets_LMs_combined.csv", "_optic_parameters.pdf", curr_file)),
                                           verbose = TRUE)
  
  # add results to tibble
  eye_L <- eye_L %>% 
    left_join(optic_parameters, by = "facet")
  
  # define viridis colours for the different parameters
  cols_to_use <- viridis(n=100, begin = 0, end = 1)
  
  eye_L <- eye_L %>% 
    arrange(facet) %>% 
    mutate(size_cols = continuous_color_ramp(eye_L$size,cols_to_use)) %>% 
    mutate(P_cols = continuous_color_ramp(eye_L$P,cols_to_use)) %>% 
    mutate(v_cols = continuous_color_ramp(eye_L$v,cols_to_use)) %>% 
    mutate(CPD_cols = continuous_color_ramp(eye_L$CPD,cols_to_use)) %>% 
    mutate(delta_phi.deg_cols = continuous_color_ramp(eye_L$delta_phi.deg,cols_to_use)) %>% 
    mutate(number_of_neighs_cols = continuous_color_ramp(eye_L$number.of.neighbours,cols_to_use))
  
  
  # # add right eye
  # eye_R <- eye_L %>% 
  #   mutate(x = -1*x,
  #          norm.x = -1*norm.x,
  #          facet = row_number() + nrow(eye_L),
  #          ID = gsub("^L_", "R_", ID))
  # 
  # # combine eyes
  # both_eyes <- rbind(eye_L,
  #                    eye_R)
  
  # save eye infos
  write_csv(as.data.frame(eye_L),
            file.path(facet_infos, gsub("_facets_LMs_combined", "_facet_infos", curr_file)))
  
  # save landmark infos
  write_csv(as.data.frame(curr_LMs),
            file.path(facet_infos, gsub("_facets_LMs_combined", "_landmarks", curr_file)))
  
  
  
  
  
  # plot eye varaibles ------------------------------------------------------  
  # # # size
  curr_variable <- "facet size"
  close3d()
  # mean size
  plot3d(both_eyes %>% 
           select(x,y,z), 
         radius = both_eyes$size,
         # size = mean(both_eyes$size)*1.2, 
         col = both_eyes$size_cols, 
         type="s",
         label = T, 
         add = F, 
         aspect = "iso")
  
  spheres3d(curr_LMs %>% select(x,y,z),
            col="blue", radius=mean(c(dist(range(eye_L$x)),dist(range(eye_L$y)),dist(range(eye_L$z))))/20)
  text3d(curr_LMs %>% select(x,y,z),
         texts = curr_LMs$ID)
  
  title3d(paste0(gsub("_rotated_facets\\.csv", "", curr_file), ": ", curr_variable),
          col = 'black', line = 5)
  
  view <- structure(c(0.800687074661255, -0.15804435312748, 0.577858686447144, 
                      0, 0.598479747772217, 0.254242211580276, -0.759724497795105, 
                      0, -0.0268459133803844, 0.954138934612274, 0.298155069351196, 
                      0, 0, 0, 0, 1), dim = c(4L, 4L))
  par3d(windowRect = c(20, 30, 1600, 1600),
        userMatrix = view)
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "1_size", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '.png', curr_output_file_name)), fmt = 'png')
  
  # draw vectors
  vec.mult <- 3
  
  plot3d(both_eyes %>% 
           select(x,y,z), 
         radius = both_eyes$size,
         # size = mean(both_eyes$size)*1.2, 
         col = both_eyes$size_cols, 
         type="s",
         label = T, 
         add = F, 
         aspect = "iso")
  
  x_start <- both_eyes %>% pull(x)
  y_start <- both_eyes %>% pull(y)
  z_start <- both_eyes %>% pull(z)
  
  x_end <- both_eyes %>% pull(x) + both_eyes %>%
    pull(norm.x)*vec.mult*(mean(both_eyes %>% pull(size)))
  y_end <- both_eyes %>% pull(y) + both_eyes %>%
    pull(norm.y)*vec.mult*(mean(both_eyes %>% pull(size)))
  z_end <- both_eyes %>% pull(z) + both_eyes %>%
    pull(norm.z)*vec.mult*(mean(both_eyes %>% pull(size)))
  
  # Combine start and end points into a single vector
  segments <- cbind(x_start, y_start, z_start, x_end, y_end, z_end)
  segments3d(x=as.vector(t(segments[,c(1,4)])),
             y=as.vector(t(segments[,c(2,5)])),
             z=as.vector(t(segments[,c(3,6)])),
             col = "purple")
  
  # for checking
  curr_eye <- both_eyes %>%
    filter(grepl("^L_", ID)) %>%
    filter(z<(-500))
  plot3d(curr_eye %>% 
           select(x,y,z),
         # radius = 7, #curr_eye$size,
         size = 15, #mean(both_eyes$size)*1.2,
         col = curr_eye$size_cols,
         # type="s",
         label = T,
         add = F,
         aspect = "iso")
  
  
  text3d(curr_eye %>%
           select(x,y,z),
         texts = both_eyes %>%
           pull(ID),
         cex=1,
         pos = 1)
  
  x_start <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(x)
  y_start <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(y)
  z_start <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(z)
  
  x_end <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(x) + both_eyes %>%
    filter(grepl("^L_", ID)) %>% pull(norm.x)*vec.mult*(mean(both_eyes %>% filter(grepl("^L_", ID)) %>% pull(size)))
  y_end <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(y) + both_eyes %>%
    filter(grepl("^L_", ID)) %>% pull(norm.y)*vec.mult*(mean(both_eyes %>% filter(grepl("^L_", ID)) %>% pull(size)))
  z_end <- both_eyes %>% filter(grepl("^L_", ID)) %>% pull(z) + both_eyes %>%
    filter(grepl("^L_", ID)) %>% pull(norm.z)*vec.mult*(mean(both_eyes %>% filter(grepl("^L_", ID)) %>% pull(size)))
  
  # Combine start and end points into a single vector
  segments <- cbind(x_start, y_start, z_start, x_end, y_end, z_end)
  segments3d(x=as.vector(t(segments[,c(1,4)])),
             y=as.vector(t(segments[,c(2,5)])),
             z=as.vector(t(segments[,c(3,6)])),
             col = "purple")
  
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  
  
  # print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "1_size", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '_vectors.png', curr_output_file_name)), fmt = 'png')
  
  # mean delta.phi.deg
  curr_variable <- "IF angle"
  close3d()
  plot3d(both_eyes %>%
           select(x,y,z),
         radius = both_eyes$size,
         # size = mean(both_eyes$size)*1.2, 
         col = eye_L$delta_phi.deg_cols,
         type="s",
         label = T,
         add = F,
         aspect = "iso")
  
  spheres3d(curr_LMs %>% select(x,y,z),
            col="blue", radius=mean(c(dist(range(eye_L$x)),dist(range(eye_L$y)),dist(range(eye_L$z))))/20)
  text3d(curr_LMs %>% select(x,y,z),
         texts = curr_LMs$ID)
  
  title3d(paste0(gsub("_rotated_facets\\.csv", "", curr_file), ": ", curr_variable),
          col = 'black', line = 5)
  
  view <- structure(c(0.800687074661255, -0.15804435312748, 0.577858686447144, 
                      0, 0.598479747772217, 0.254242211580276, -0.759724497795105, 
                      0, -0.0268459133803844, 0.954138934612274, 0.298155069351196, 
                      0, 0, 0, 0, 1), dim = c(4L, 4L))
  par3d(windowRect = c(20, 30, 1600, 1600),
        userMatrix = view)
  
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  # print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "2_IF_angle", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '.png', curr_output_file_name)), fmt = 'png')
  
  
  # mean P
  curr_variable <- "eye parameter (P)"
  close3d()
  
  plot3d(both_eyes %>%
           select(x,y,z),
         radius = both_eyes$size,
         # size = mean(both_eyes$size)*1.2, 
         col = eye_L$P_cols,
         type="s",
         label = T,
         aspect = "iso")
  
  spheres3d(curr_LMs %>% select(x,y,z),
            col="blue", radius=mean(c(dist(range(eye_L$x)),dist(range(eye_L$y)),dist(range(eye_L$z))))/20)
  text3d(curr_LMs %>% select(x,y,z),
         texts = curr_LMs$ID)
  
  title3d(paste0(gsub("_rotated_facets\\.csv", "", curr_file), ": ", curr_variable),
          col = 'black', line = 5)
  
  view <- structure(c(0.800687074661255, -0.15804435312748, 0.577858686447144, 
                      0, 0.598479747772217, 0.254242211580276, -0.759724497795105, 
                      0, -0.0268459133803844, 0.954138934612274, 0.298155069351196, 
                      0, 0, 0, 0, 1), dim = c(4L, 4L))
  par3d(windowRect = c(20, 30, 1600, 1600),
        userMatrix = view)
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  # print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "3_P", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '.png', curr_output_file_name)), fmt = 'png')
  
  
  # mean CPD
  curr_variable <- "cycles per degree (CPD)"
  close3d()
  plot3d(both_eyes %>%
           select(x,y,z),
         radius = both_eyes$size,
         # size = mean(both_eyes$size)*1.2, 
         col = eye_L$CPD_cols,
         type="s",
         label = T,
         aspect = "iso")
  
  spheres3d(curr_LMs %>% select(x,y,z),
            col="blue", radius=mean(c(dist(range(eye_L$x)),dist(range(eye_L$y)),dist(range(eye_L$z))))/20)
  text3d(curr_LMs %>% select(x,y,z),
         texts = curr_LMs$ID)
  
  title3d(paste0(gsub("_rotated_facets\\.csv", "", curr_file), ": ", curr_variable),
          col = 'black', line = 5)
  
  view <- structure(c(0.800687074661255, -0.15804435312748, 0.577858686447144, 
                      0, 0.598479747772217, 0.254242211580276, -0.759724497795105, 
                      0, -0.0268459133803844, 0.954138934612274, 0.298155069351196, 
                      0, 0, 0, 0, 1), dim = c(4L, 4L))
  par3d(windowRect = c(20, 30, 1600, 1600),
        userMatrix = view)
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  # print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "4_CPD", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '.png', curr_output_file_name)), fmt = 'png')
  
  # of neighbours
  curr_variable <- "# of neighbours"
  close3d()
  plot3d(eye_L %>%
           select(x,y,z),
         radius = eye_L$size,
         col = eye_L$number_of_neighs_cols,
         type="s",
         label = T,
         aspect = "iso")
  
  title3d(paste0(gsub("_rotated_facets\\.csv", "", curr_file), ": ", curr_variable),
          col = 'black', line = 5)
  
  view <- structure(c(0.800687074661255, -0.15804435312748, 0.577858686447144, 
                      0, 0.598479747772217, 0.254242211580276, -0.759724497795105, 
                      0, -0.0268459133803844, 0.954138934612274, 0.298155069351196, 
                      0, 0, 0, 0, 1), dim = c(4L, 4L))
  par3d(windowRect = c(20, 30, 1600, 1600),
        userMatrix = view)
  
  bbox3d(alpha = 0.0, xlab="NULL")
  bg3d(sphere = TRUE, color = "transparent")
  
  # print("Saving rotated data as png...")
  curr_output_file_name <- gsub("rotated_facets", "0_no_of_neighbors", curr_file)
  # rgl.postscript(file.path(facet_infos, paste0(curr_output_file_name, '.pdf')), fmt = 'svg')
  rgl.snapshot(file.path(facet_infos, gsub('.csv', '.png', curr_output_file_name)), fmt = 'png')
  
  
  # 
  # # draw vectors
  # vec.mult <- mean(eye_L$size)*10
  # for(curr_facet in 1:nrow(eye_L)){ # nrow(eye_L)
  #   normal_vectors_df_subset <- eye_L %>%
  #     filter(facet == curr_facet) %>%
  #     select(norm.x, norm.y, norm.z)
  #   curr_facet_coordinates <- eye_L %>% filter(facet==curr_facet) %>% select(x,y,z)
  # 
  #   # find mean point of normalized normal vector ends
  #   norm.y <- normal_vectors_df_subset$norm.y
  #   norm.z <- normal_vectors_df_subset$norm.z
  # 
  #   lines3d(x = c(curr_facet_coordinates %>% pull(x), curr_facet_coordinates %>% pull(x) + norm.x*vec.mult),
  #           y = c(curr_facet_coordinates %>% pull(y), curr_facet_coordinates %>% pull(y) + norm.y*vec.mult),
  #           z = c(curr_facet_coordinates %>% pull(z), curr_facet_coordinates %>% pull(z) + norm.z*vec.mult),
  #           col = eye_L %>%
  #             filter(facet == curr_facet) %>%
  #             select(delta_phi.deg_cols ))
  # }
  print("*******************************")
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
                                 facet_estimate = 14,
                                 cores = 18,
                                 plot_results = FALSE,
                                 verbose = FALSE){
  
  # # testing
  # df = LMs_facets_combined
  # plot_file = file.path(facet_infos_folder,
  #                       gsub("_surface.stl", "_neighbour_and_size_data.pdf", basename(file_name)))
  # facet_estimate = 14
  # cores = 18
  # plot_results = TRUE
  # verbose = TRUE
  
  # get data of facets
  curr_facets <- df %>% 
    filter(type == "facet")
  
  # get data of LMs
  curr_LMs <- df %>% 
    filter(type == "LM")
  
  # # get data for left eye
  # eye_L <- curr_facets %>% 
  #   mutate(facet = row_number(),
  #          ID = paste0("L_", ID)) %>%
  #   select(c(ncol(.), 1:(ncol(.)-1)))
  
  # # do the next calculations only for left eye and mirror results to right eye
  
  df_neighbours <- find_neighbours(df = curr_facets,
                                   facet_size = facet_estimate,
                                   neighbour_threshold = 1.5,
                                   cores = cores,
                                   plot_results = plot_results,
                                   plot_file = plot_file,
                                   verbose = verbose)
  
  
  
  # add results to tibble
  df_w_neighbours <- curr_facets %>% 
    left_join(df_neighbours, 
              by = "ID")
  
  
  # calculate facet normals according to their neighbors
  df_normals <- get_facet_normals(df = df_w_neighbours,
                                  cores = cores,
                                  plot_file = gsub("_neighbour_and_size_data", 
                                                   "_normal_data", 
                                                   plot_file),
                                  verbose = TRUE)
  
  
  # add info to df
  df_w_normals <- df_w_neighbours %>% 
    left_join(df_normals, 
              by="ID")
  
  
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
  
  return(df_w_plot_cols)
}