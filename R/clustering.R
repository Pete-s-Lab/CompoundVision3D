#' Find threshold for agglomerative clustering
#'
#' Very roughly find clusters around facet peaks.
#'
#' find good threshold - choose x,y, x,z or other dimension combinations to get reasonable plots
#' @export
find_threshold <- function(df,
                           column1 = NULL,
                           column2 = NULL,
                           height_column,
                           min_threshold = 1,
                           max_treshold = 2.5, 
                           trials = 9,
                           plot_file = NULL){
  
  # if undefined, take the columns with the highest ranges for plotting
  if(is.null(column1)){
    column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
    names(column_ranges) <- c("x","y","z")
    column1 <- names(sort(column_ranges))[2]
  }
  if(is.null(column2)){
    column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
    names(column_ranges) <- c("x","y","z")
    column2 <- names(sort(column_ranges))[3]
  }
  
  if(trials > 4){
    number_y = 3
  } else if(trials > 2) {
    number_y = 2
  } else{
    number_y = 1
  }
  number_x = ceiling(trials / number_y)
  
  par(mfrow=c(number_x,number_y))
  for(curr_thrsh in round(seq(min_threshold, max_treshold, length.out = trials),2)){
    df_tmp <- df %>% 
      filter(!!as.symbol(height_column) >= curr_thrsh) 
    plot(df_tmp %>%
           select(!!as.symbol(column1),
                  !!as.symbol(column2)),
         col = "black",
         pch=16,
         cex = .01,
         main=paste0(round(curr_thrsh,2), "; n = ", nrow(df_tmp)))
  }
  par(mfrow=c(1,1))
  
  if(!is.null(plot_file)){
    pdf(file = plot_file, # , today()
        onefile = TRUE, width = (21.0-4)/1.5, height = (((21.0-4)))*3/4)
    par(mfrow=c(3,3))
    for(curr_thrsh in round(seq(min_threshold, max_treshold, length.out = trials),2)){
      df_tmp <- local_heights %>% 
        filter(!!as.symbol(height_column) >= curr_thrsh) 
      plot(df_tmp %>%
             select(!!as.symbol(column1),
                    !!as.symbol(column2)),
           col = "black",
           pch=16,
           cex = .1,
           main=paste0(round(curr_thrsh,2), "; n = ", nrow(df_tmp)),
           asp = 1)
    }
    dev.off()
    par(mfrow=c(1,1))
  }
}



#' Rough clustering towards Facet Peaks
#'
#' Very roughly find clusters around facet peaks.
#'
#' @param df A tibble containing triangle center coordinates of vertices that 
#' lie above threshold in columns `x, y, z` and local height in `local_height`. 
#' Typically, this is the resulting tibble of the `get_local_height()` function.
#' @param local_height_threshold Threshold to filter out vertices: Vertices with
#' `local_height >= local_height_threshold * mean(local_height)` will be kept. 
#' Default: `2.5`.
#' @param clust_melt_rad A numerical value of search radius for agglomerative 
#' clustering.
#' @param iterations A numerical value of how many clustering iterations should 
#' be run. Default: `1`.
#' @param cores A numerical value of how many cores to use. Default: `1`.
#' @return Tibble `df` with ideally one vertex per facet, located on the facet 
#' peak.
#'
#' @export
#' @examples
#' xxx: add example
#' 
find_facets_rough <- function(df,
                              local_height_threshold = 2.5,
                              height_column,
                              clust_melt_rad, 
                              # iterations = 1,
                              cores = 1){
  
  # # testing
  # local_height_threshold = 2
  # df = local_heights
  # height_column = "local_height_log"
  # clust_melt_rad = 1.0*round(1/8*curr_search_diam,2)
  # cores = 12
  
  # load multi-core package
  require(doParallel)
  
  start_time <- Sys.time()
  
  # dplyr NULLs
  ID <- x <- y <- z <- i <- local_height <- NULL
  
  # remove all coordinates with a local height of less than 2.5x mean local height
  df <- df %>% 
    filter(!!as.symbol(height_column) >= local_height_threshold) 
  
  
  # # store filtered df for final results
  # df_final <- df
  
  # select relevant columns for analyses within this function
  df <- df %>% 
    select(x, y, z)
  
  # # plot filtered tibble in 'SEM colours'
  # plot3d(df,
  #        aspect = "iso")
  
  message("Getting agglomerative clusters using ", nrow(df), " vertices to start with ...")
  
  registerDoParallel(cores)
  # for(k in 1:iterations){
  # print(paste("iteration #:", k))
  # define number of iterations of preliminary agglomerative clustering
  local.clust.verts.means <- foreach(i = 1:nrow(df),
                                     .combine=rbind, .packages=c('dplyr')) %dopar% { # .packages=c('dplyr', 'geometry')
                                       
                                       # get last columns of df because those are the ones containing the results of the last iteration
                                       curr.df <- df[,(ncol(df)-2):ncol(df)]
                                       colnames(curr.df) <- c("x", "y", "z")
                                       
                                       local.clust.verts <- curr.df %>%
                                         filter(x>(df$x[i]-clust_melt_rad) & x<(df$x[i]+clust_melt_rad) &
                                                  y>(df$y[i]-clust_melt_rad) & y<(df$y[i]+clust_melt_rad) &
                                                  z>(df$z[i]-clust_melt_rad) & z<(df$z[i]+clust_melt_rad))
                                       # plot(local.clust.verts[,2:3])
                                       tmp <- local.clust.verts %>% # local.clust.verts.mean
                                         summarize(median_x = median(x),
                                                   median_y = median(y),
                                                   median_z = median(z))
                                       # points(local.clust.verts.mean[,2:3], pch=16,col="red")
                                     }
  
  # colnames(local.clust.verts.means) <- as.vector(outer(c("local_meds_x", "local_meds_y", "local_meds_z"), 
  #                                                      k, 
  #                                                      paste, sep="_"))
  # df_final <- as_tibble(cbind(df_final, local.clust.verts.means))
  
  # dist.matr <- dist(local.clust.verts.means)
  # dist.matr.tbl <- melt(as.matrix((dist.matr))) %>% as_tibble() %>% filter(Var1 < Var2) %>% 
  #   filter(value <= 50)
  # png(paste0(stl.folder, "/hist_below_50_",k,".png"))
  # hist(dist.matr.tbl$value)
  # dev.off()
  
  # }
  stopImplicitCluster()
  
  
  print("done!")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  return(local.clust.verts.means)
}




#' Fine clustering towards Facet Peaks
#'
#' Find facet peaks.
#'
#' @param df A tibble containing triangle center coordinates of vertices that 
#' lie above threshold in columns `x, y, z` . 
#' Typically, this is the resulting tibble of the `find_facets_rough()` 
#' function.
#' @param cols_to_use A numerical value of which columns store x,y,z coordinates 
#' of rough clustering.
#' @return Tibble `df` with ideally one vertex per facet, located on the facet 
#' peak.
#'
#' @export
#' @examples
#' xxx: add example
#' 

find_facet_canidates <- function(df,
                                 cols_to_use = 1:3,
                                 h_min = NULL,
                                 h_max = NULL,
                                 h_final = NULL,
                                 column1 = NULL,
                                 column2 = NULL,
                                 n_steps = 100,
                                 trials = 9,
                                 plot_file = NULL,
                                 verbose = FALSE){
  
  # # testing
  # df = facet_candidates
  # cols_to_use = 1:3
  # h_min = 20
  # h_max = 8
  # h_final = 8
  # column1 = NULL
  # column2 = NULL
  # column1 = "x"
  # column2 = "y"
  # n_steps = 100
  # plot_file = NULL
  # verbose = TRUE
  
  # set up function to plot
  cadidate_plotting <- function(){
    # plot dendrogram
    x_total <- length(hc1$order)
    y_total <- max(hc1$height)
    
    # set x and y limits of dendrogram plot to 1/4 and 1/5th, respecitivels
    xlim_frac <- c(0, x_total / 4)
    ylim_frac <- c(0, y_total / 5)
    
    # if there are fewer branches than 90, set to 180 or to total number
    if((x_total / 4) < 90) xlim_frac[2] <- 90
    if(x_total < 90) xlim_frac[2] <- x_total
    
    
    # Set plot layout
    layout(mat = matrix(c(1, 2, 3, 4), 
                        nrow = 4, 
                        ncol = 1),
           heights = c(2,1,1,3))# ,    # Heights of the two rows
    
    # plot dendrogram
    op <- par(mar = c(1, 4, 3, 2) + 0.1)  # bottom, left, top, right par(mar = c(1, 4, 3, 2) + 0.1)
    plot(
      as.dendrogram(hc1),
      cex = 0.1,
      xlab = "",
      main = paste0("Zoom of hierarchical clustering dendrogram (cur. cf.: ", 
                    as.numeric(curr_cutoff), 
                    "; final cf.: ", as.numeric(h_final$x),
                    ")"),
      sub = "",
      xlim = xlim_frac,
      ylim = ylim_frac,
      leaflab = "none"
    )
    par(op)
    
    abline(a=h_max, b=0, col="blue", lty=2)
    abline(a=h_min, b=0, col="blue", lty=2)
    abline(a=as.numeric(curr_cutoff), b=0, col="blue", lty=1)
    abline(a=as.numeric(h_final), b=0, col="red", lty=1)
    
    
    # curve
    plot(ommatidia.no.df$h, ommatidia.no.df$ommatidia.no,
         xlab = "Cutoff value",
         ylab = "facet number") # , ylim = c(0,max(ommatidia.no.df$ommatidia.no))
    lines(x=rep(h_min, 2), y = c(min(ommatidia.no.df$ommatidia.no), max(ommatidia.no.df$ommatidia.no)),
          col = "blue", lty=2)
    lines(x=rep(h_max, 2), y = c(min(ommatidia.no.df$ommatidia.no), max(ommatidia.no.df$ommatidia.no)),
          col = "blue", lty=2)
    lines(x=rep(as.numeric(curr_cutoff$x), 2), y = c(min(ommatidia.no.df$ommatidia.no), max(ommatidia.no.df$ommatidia.no)),
          col = "blue")
    lines(x=rep(as.numeric(h_final$x), 2), y = c(min(ommatidia.no.df$ommatidia.no), max(ommatidia.no.df$ommatidia.no)),
          col = "red")
    
    # differences
    plot(ommatidia.no.df$h, ommatidia.no.df$ommatidia.no.diff, type="l",
         xlab = "Cutoff value",
         ylab = "Delta facet number")
    lines(x=rep(h_min, 2), y = c(min(ommatidia.no.df$ommatidia.no.diff), max(ommatidia.no.df$ommatidia.no.diff)),
          col = "blue", lty=2)
    lines(x=rep(h_max, 2), y = c(min(ommatidia.no.df$ommatidia.no.diff), max(ommatidia.no.df$ommatidia.no.diff)),
          col = "blue", lty=2)
    lines(x=rep(as.numeric(curr_cutoff$x), 2), y = c(min(ommatidia.no.df$ommatidia.no.diff), max(ommatidia.no.df$ommatidia.no.diff)),
          col = "blue")
    lines(x=rep(as.numeric(h_final$x), 2), y = c(min(ommatidia.no.df$ommatidia.no.diff), max(ommatidia.no.df$ommatidia.no.diff)),
          col = "red")
    abline(a=0, b=0, col="red", lty=2)
    
    
    plot(df.fin.clean %>% 
           select(!!column1, !!column2),
         asp = 1,
         pch = 16,
         cex = .5)
    
    # # Histogram
    # hist.plot <- hist(dist.clusters.tbl$value,
    #                   main = "Histogram of distances",
    #                   xlab = paste0("Distance"))
    # break_size <- hist.plot$breaks[2]
    # hist_x <- hist.plot$breaks[which(hist.plot$counts == max(hist.plot$counts))[1]] + break_size/2
    # lines(x = rep(hist_x, 2), y=c(0, (max(hist.plot$counts) + 0.05 * max(hist.plot$counts))),
    #       col = "blue")
    # par(mfrow=c(1,1))
  }
  
  # if undefined, take the columns with the highest ranges for plotting
  if(is.null(column1)){
    column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
    names(column_ranges) <- c("x","y","z")
    column1 <- names(sort(column_ranges))[2]
  }
  if(is.null(column2)){
    column_ranges <- c(diff(range(df$x)), diff(range(df$y)), diff(range(df$z)))
    names(column_ranges) <- c("x","y","z")
    column2 <- names(sort(column_ranges))[3]
  }
  
  start_time <- Sys.time()
  # dplyr NULLs
  ID <- x <- y <- z <- cluster <- Var1 <- Var2 <- '.' <- NULL
  
  # Agglomerative clustering
  # Dissimilarity matrix
  if(verbose == TRUE){
    cat("Calculating distance matrix...", "\n")
  }
  d <- dist(df[,cols_to_use], method = "euclidean")
  hc1 <- hclust(d, method = "complete" )
  
  if(is.null(h_min) | is.null(h_max)){
    # Hierarchical clustering using Complete Linkage
    
    if(verbose == TRUE){
      cat("Calculating and plotting dendrogram of hierarchichal clustering...", "\n")
    }
    
    # Plot the obtained dendrogram
    cat("select minimum and maximum cut-off points on y axis for first trial.", "\n")
    cat("Only select two points, the script will continue automatically.", "\n")
    
    h.cutoff.df <- locator(type = "n", n=2)
    h_min = h.cutoff.df$y[1]
    h_max = h.cutoff.df$y[2]
    
    if(verbose == TRUE) cat("Cutoff values:", paste0(round(h_min,3), "; ", round(h_max, 3),"."), "\n")
  } else {
    if(verbose == TRUE) cat("Cut-offs defined as:", round(h_min,3), "and", round(h_max, 3),".", "\n")
  }
  
  if(verbose == TRUE) cat(paste0("Min.: ", round(h_min,3), "; max.: ", round(h_max, 3)), "\n")
  
  # n_steps = 200
  names = c("h", "ommatidia.no")
  
  
  if(is.null(h_final)){
    cat("Select cut-off point on x-axis.", "\n")
    h_final <- locator(type = "n", n=1)
    if(verbose==TRUE) cat(paste0("Final cut-off chosen: ", round(h_final$x[length(h_final$x)], 2)), "\n")
  } else if(is.numeric(h_final)){
    h_final <- data.frame(x=h_final)
    if(verbose==TRUE) cat(paste0("Final cut-off value defined as: ", round(h_final,3),"."), "\n")
  }
  
  if(verbose == TRUE) cat("Finding clusters for", n_steps, "points between cut-off values...", "\n")
  
  ommatidia.no.df = as_tibble(setNames(data.frame(matrix(nrow = 0, 
                                                         ncol = length(names))), 
                                       names))
  
  for(h in seq(h_min, h_max, length.out = n_steps)){
    clusters.indeces <- cutree(hc1, h = h)
    ommatidia.no <- length(unique(clusters.indeces))
    ommatidia.no.df <- rbind(ommatidia.no.df, cbind(h, ommatidia.no))
    # show.progress(h, h_max)
  }
  
  ommatidia.no.df <- ommatidia.no.df %>%
    mutate(ommatidia.no.diff = ommatidia.no - lag(ommatidia.no, default = ommatidia.no[2]))
  
  
  
  # define range in which to look for cutoffs
  cutoff_range <- seq(h_max,h_min,length.out=trials)
  cutoff_range <- cutoff_range[cutoff_range >= 1]
  
  # plot ranges
  # save to file
  if(!is.null(plot_file)){
    # curr_plot_file <- gsub("_facet_candidates\\.pdf", 
    #                        paste0("_facet_candidates_", 
    #                               as.numeric(curr_cutoff),
    #                               ".pdf"), 
    #                        plot_file)
    if(verbose == TRUE) cat("Saving plots as ", plot_file, "\n")
    
    # PDF plots
    pdf(plot_file, # , today()
        onefile = TRUE, 
        paper = "a4", 
        height = 14)
  }
  
  curr_cutoff = cutoff_range[1]
  ctr = 0
  for(curr_cutoff in cutoff_range){
    ctr = ctr+1
    if(verbose==TRUE) cat("Current cutoff:", curr_cutoff, 
                          "(", ctr, "/", length(cutoff_range), ")\n")
    
    curr_cutoff <- data.frame(x=curr_cutoff)
    
    
    # save a vector (clusters.fin) that stores to which cluster each coordinate belongs
    clusters.fin <- cutree(hc1, h = curr_cutoff$x[length(curr_cutoff$x)])
    
    # add cluster-memberships (clusters.fin) to the cluster coordinates (df)
    # reduce clusters to their mean coordinates
    df.fin <- df %>% 
      select(all_of(cols_to_use)) %>% 
      mutate(cluster = clusters.fin) %>% 
      group_by(cluster) %>% 
      dplyr::summarize(x = median(!!as.symbol(colnames(df)[cols_to_use[1]])),
                       y = median(!!as.symbol(colnames(df)[cols_to_use[2]])),
                       z = median(!!as.symbol(colnames(df)[cols_to_use[3]]))) %>% 
      select(-cluster)
    
    # # create distance matrix of all clusters to each other
    # dist.clusters <- dist(df.fin)
    # 
    # # melt the distance matrix into a three-column tibble
    # dist.clusters.tbl <- reshape2::melt(as.matrix((dist.clusters))) %>% as_tibble() %>% filter(Var1 < Var2)
    # 
    # clust.dist.med <- median(dist.clusters.tbl$value)
    
    # filter the clusters that remain
    df.fin.clean <- df.fin %>% 
      mutate(ID = 1:nrow(.))
    
    
    
    # # device plotting
    cadidate_plotting()
  }
  
  if(!is.null(plot_file)){
    dev.off()
  }
  
  # reset to h_final
  if(verbose==TRUE) cat("Resetting to:", as.numeric(h_final), "\n")
  
  curr_cutoff <- h_final
  
  
  # save a vector (clusters.fin) that stores to which cluster each coordinate belongs
  clusters.fin <- cutree(hc1, h = curr_cutoff$x[length(curr_cutoff$x)])
  
  # add cluster-memberships (clusters.fin) to the cluster coordinates (df)
  # reduce clusters to their mean coordinates
  df.fin <- df %>% 
    select(all_of(cols_to_use)) %>% 
    mutate(cluster = clusters.fin) %>% 
    group_by(cluster) %>% 
    dplyr::summarize(x = median(!!as.symbol(colnames(df)[cols_to_use[1]])),
                     y = median(!!as.symbol(colnames(df)[cols_to_use[2]])),
                     z = median(!!as.symbol(colnames(df)[cols_to_use[3]]))) %>% 
    select(-cluster)
  
  # # create distance matrix of all clusters to each other
  # dist.clusters <- dist(df.fin)
  # 
  # # melt the distance matrix into a three-column tibble
  # dist.clusters.tbl <- reshape2::melt(as.matrix((dist.clusters))) %>% as_tibble() %>% filter(Var1 < Var2)
  # 
  # clust.dist.med <- median(dist.clusters.tbl$value)
  
  # filter the clusters that remain
  df.fin.clean <- df.fin %>% 
    mutate(ID = 1:nrow(.))
  
  if(verbose == TRUE){
    cat(paste0("Found ", nrow(df.fin.clean), " facet center candiates. Check 3D plot device."), "\n")
    
    end_time <- Sys.time()
    cat("Time taken for manual and anutomatic process:", end_time - start_time, "\n")
  }
  
  # plot eye in 3D to get overview
  if(verbose == TRUE){
    cat("Plotting 'SEM'-coloured eye in RGL 3D window. Check it out to get overview.", "\n")
  }
  plot3d(df %>%
           select(all_of(cols_to_use)),
         # col = local_heights$local_height_col_log,
         aspect = "iso")
  
  # plot the cluster centers
  spheres3d(df.fin.clean %>%
              select(x,y,z),
            col="red", radius=5, alpha = 1)
  
  return(df.fin.clean %>% 
           mutate(cutoff_min = round(h_min,5),
                  cutoff_max = round(h_max, 5),
                  cutoff_fin = round(as.numeric(h_final$x), 5)))
}




#' Agglomerative clustering
#' 
#' xxx: add description
#' 
#' @param df A tibble containing 3D point coordinates of vertices that 
#' #' lie above threshold. Should contain the columns `x, y, z`. 
#' #' Typically, this is the resulting tibble of the `get_local_height()` function.
#' @param clust_melt_rad A numerical value of search radius for agglomerative
#' clustering.
#' @param iterations A numerical value of how many clustering iterations should
#' be run. Default: `1`.
#' @param cores Number of CPU cores to use. Default = `2`.
#'
#' @export
#' @examples
#' xxx: add example
#' 
agglomerative_clustering <- function(df,
                                     clust_melt_rad,
                                     iterations,
                                     cores = 2){
  # # testing
  # df = rough_peaks %>%
  # slice(1:1000)
  # clust_melt_rad = (curr_search_diam/128*counter/2)
  # iterations = 1
  # cores=12
  
  # dplyr NULLs
  ID <- x <- y <- z <- i <- local_height <- NULL
  
  df <- df %>% 
    select(x, y, z)
  
  registerDoParallel(cores)
  
  k=1
  i=1
  for(k in 1:iterations){
    print(paste0("iteration #: ",k))
    if(k==1){
      df_final <- df
      print(paste0("original number of vertices: ", nrow(df_final)))
    }
    
    curr_clust_melt_rad <- clust_melt_rad # clust_melt_rad + ((0) * 0.1) * clust_melt_rad # k-1
    print(paste0("Current melting radius: ", round(curr_clust_melt_rad, 3)))
    
    
    # define number of iterations of preliminary agglomerative clustering
    local.clust.verts.means <- foreach(i = 1:nrow(df_final),
                                       .combine=rbind, .packages=c('dplyr')) %dopar% { # .packages=c('dplyr', 'geometry')
                                         
                                         local.clust.verts <- df_final %>%
                                           filter(x>(df_final$x[i]-curr_clust_melt_rad) & x<(df_final$x[i]+curr_clust_melt_rad) &
                                                    y>(df_final$y[i]-curr_clust_melt_rad) & y<(df_final$y[i]+curr_clust_melt_rad) &
                                                    z>(df_final$z[i]-curr_clust_melt_rad) & z<(df_final$z[i]+curr_clust_melt_rad))
                                         # plot(local.clust.verts[,2:3])
                                         tmp <- local.clust.verts %>% # local.clust.verts.mean
                                           dplyr::summarize(median_x = median(x),
                                                            median_y = median(y),
                                                            median_z = median(z))
                                         # points(local.clust.verts.mean[,2:3], pch=16,col="red")
                                       }
    
    colnames(local.clust.verts.means) <- c("x", "y", "z")
    
    # colnames(local.clust.verts.means) <- as.vector(outer(c("local_meds_x", "local_meds_y", "local_meds_z"),
    #                                                      k,
    #                                                      paste, sep="_"))
    
    df_final <- 
      distinct(local.clust.verts.means %>% 
                 round())
    
    
    
    # df_final <- local.clust.verts.means
    
    
    print(paste0("reduced number of vertices: ", nrow(df_final)))
    
    # dist.matr <- dist(local.clust.verts.means)
    # dist.matr.tbl <- melt(as.matrix((dist.matr))) %>% as_tibble() %>% filter(Var1 < Var2) %>% 
    #   filter(value <= 50)
    # png(paste0(stl.folder, "/hist_below_50_",k,".png"))
    # hist(dist.matr.tbl$value)
    # dev.off()
    
  }
  stopImplicitCluster()
  
  # plot3d(df, size = 3)
  # points3d(df_final, col="red", size = 5)
  
  return(df_final)
}