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
                                    verbose = FALSE){
  
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
                           .combine=rbind, .packages=c('dplyr', 'geometry')) %dopar% {
                             
                             curr.facet.x.y.z <- df %>%
                               dplyr::slice(i) %>%
                               select(x, y, z) %>%
                               as.numeric()
                             
                             search_diam_local_height <- round(1/8*search_diam,2)
                             
                             curr.facets.df <- df %>%
                               dplyr::filter(x  >= curr.facet.x.y.z[1] - search_diam_local_height &
                                               y  >= curr.facet.x.y.z[2] - search_diam_local_height &
                                               z  >= curr.facet.x.y.z[3] - search_diam_local_height &
                                               x  <= curr.facet.x.y.z[1] + search_diam_local_height &
                                               y  <= curr.facet.x.y.z[2] + search_diam_local_height &
                                               z  <= curr.facet.x.y.z[3] + search_diam_local_height )
                             # print(nrow(curr.facets.df))
                             
                             # calculate current average facet normal
                             curr.facets.av.normal <- c(mean(curr.facets.df$norm.x), 
                                                        mean(curr.facets.df$norm.y), 
                                                        mean(curr.facets.df$norm.z))
                             
                             # calculate current facet center
                             curr.facets.center <- c(mean(curr.facets.df$x), 
                                                     mean(curr.facets.df$y), 
                                                     mean(curr.facets.df$z))
                             
                             
                             # create unit normal vector of plane normal
                             curr.facets.av.normal.normailzed <- normalize_vector(curr.facets.av.normal)
                             
                             # find vector of current point to arbitrary other point on plane (here: plane center)
                             vector_point_facet_center <- curr.facet.x.y.z-curr.facets.center
                             
                             
                             curr_local_height <- dot(vector_point_facet_center,curr.facets.av.normal.normailzed, 
                                                      d=TRUE)
                             
                             # plot3d(curr.facets.df[2:4], aspect = "iso")
                             # par3d("windowRect"= c(2300,200,3400,1000))
                             # lines3d(x=c(curr.facets.center[1], curr.facets.center[1]+curr.facets.av.normal.normailzed[1]),
                             #         y=c(curr.facets.center[2], curr.facets.center[2]+curr.facets.av.normal.normailzed[2]),
                             #         z=c(curr.facets.center[3], curr.facets.center[3]+curr.facets.av.normal.normailzed[3]))
                             
                             # df$local_height[n] <- curr_local_height
                             
                             tmp <- curr_local_height
                           }
  
  # add distances to local planes to df tibble
  df$local_height <- as.numeric(local_heights)
  # save(df, file = "./data/AV00003_df.Rdata")
  
  stopImplicitCluster()
  
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
  
  # limet number of rows for one analsis to get reasonable calulation times
  max_rows = 2500
  starts <- seq(1,nrow(df),max_rows)
  
  if(verbose == TRUE){
    if(length(starts) > 1){
      cat("Splitting analysis into ", length(starts), " parts.\n")
    }
  }
  
  
  df_normalized_raw <- tibble(n=numeric(),
                              local_heights_quantiles_normalized=numeric())
  
  if(verbose == TRUE){
    cat("Starting analyses on cluster...\n")
    start_time <- Sys.time()
  }
  
  registerDoParallel(cores)
  s=1
  for(k in 1:length(starts)){
    s=starts[k]
    
    e = (s + max_rows - 1)
    if(e > nrow(df)) e = nrow(df)
    # print(paste0(s, " to ", e))
    i=s
    curr_df_normalized_raw <- foreach(i = s:e,# nrow(df)
                                      .combine=rbind,
                                      .packages=c('dplyr', 'forceR')) %dopar% {
                                        
                                        curr_vertex <- df %>%
                                          slice(i) %>%
                                          select(x,y,z)
                                        
                                        
                                        # spheres3d(curr_vertex,
                                        #           col="blue",
                                        #           radius = 10)
                                        
                                        curr_facetsized_ROI <- df %>%
                                          filter(x >= curr_vertex$x - normalize_diam,
                                                 x <= curr_vertex$x + normalize_diam,
                                                 y >= curr_vertex$y - normalize_diam,
                                                 y <= curr_vertex$y + normalize_diam,
                                                 z >= curr_vertex$z - normalize_diam,
                                                 z <= curr_vertex$z + normalize_diam) %>%
                                          select(n,x,y,z,one_of(column_to_normalize)) %>%
                                          rename(col_to_norm = 5) # rename fith column
                                        
                                        # points3d(curr_facetsized_ROI %>%
                                        #            select(x,y,z),
                                        #          col="red", size=11)
                                        
                                        # get local heights values and filter out outliers
                                        curr_local_heights <- curr_facetsized_ROI %>%
                                          select(n, col_to_norm)
                                        
                                        if(nrow(curr_local_heights)>1){
                                          # set outliers as quantile values
                                          curr_Q <- quantile(curr_local_heights$col_to_norm, probs=c(.10, .90), na.rm = FALSE)
                                          curr_Q_min <- curr_Q[1]
                                          curr_Q_max <- curr_Q[2]
                                          
                                          curr_local_heights <- curr_local_heights %>%
                                            mutate(local_heights_quantiles = case_when(col_to_norm < curr_Q_min ~ curr_Q_min,
                                                                                       col_to_norm > curr_Q_max ~ curr_Q_max,
                                                                                       TRUE ~ col_to_norm),
                                                   local_heights_quantiles_normalized = rescale_to_range(local_heights_quantiles, 0, 1))
                                          
                                          # hist(curr_local_heights$col_to_norm, breaks = 7)
                                          # hist(curr_local_heights$local_heights_quantiles, breaks = 7)
                                          # hist(curr_local_heights$local_heights_quantiles_normalized, breaks = 7)
                                          
                                          tmp <- curr_local_heights %>%
                                            select(n,local_heights_quantiles_normalized) # %>%
                                          # mutate(iteration = i)# %>%
                                          # rename_at(vars(all_of(col.from)), ~col.to)
                                        }
                                      }
    # if(verbose == TRUE){
    #   print("merging data...")
    # }
    
    df_normalized_raw <- rbind(df_normalized_raw, curr_df_normalized_raw)
    
    if(verbose == TRUE){
      cat(round(k*100/length(starts),2), "%\n")
      # print("*******")
    }
  }
  
  # terminate cluster
  stopImplicitCluster()
  
  if(verbose == TRUE){
    cat("Multi-threaded analysis finished.\n")
    end_time <- Sys.time()
    cat(end_time - start_time, "\n")
  }
  
  if(verbose == TRUE){
    cat("Summarizing ", nrow(df_normalized_raw), " results...\n")
  }
  
  # calculate means per facet - i.e. row = n
  df_normalized_summarized <- df_normalized_raw %>%
    group_by(n) %>%
    summarise(local_height_norm = mean(local_heights_quantiles_normalized ))
  
  df_fin <- df %>%
    left_join(df_normalized_summarized, by = "n") %>%
    select(-n)
  
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
             col = df_fin %>%
               pull(local_height_log_norm_col),
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
                                   facet_estimate = 14,
                                   cores,
                                   plot_results = FALSE,
                                   verbose = FALSE){
  # # testing
  # df = facet_candidates
  # local_heights = local_heights
  # landmarks_file = landmarks_file
  # crop_log_file = crop_log_file
  # facet_estimate = 14
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
    
    ROIs <- c("head", "eye1", "eye2")
    
    
    curr_ROI_coordinates <- get_ROI_coordinates(ROIs,
                                                crop_log_data)
    
    
    # add eye ROI coordinates to curr. eye ROI
    local_heights_translated <- translate_ROIs(df = local_heights,
                                               ROI_coordinates = curr_ROI_coordinates,
                                               eye = as.numeric(curr_eye))
    
    
    
    facet_positions_translated <- translate_ROIs(df = facet_positions_new,
                                                 ROI_coordinates = curr_ROI_coordinates,
                                                 eye = as.numeric(curr_eye))
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
