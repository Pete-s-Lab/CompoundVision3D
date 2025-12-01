#' Read multiple csv files and store their file name in column
#'
#' @param file_list A vector containing `character` strings of csv file names.
#' @param show_col_types A `logic` value to define if column types should be 
#' shown while importing. `Default: FALSE`.
#' @param axis A character string defining the global axis to align to. Must be 
#' `x`, `y`, or `z`.
#' @return Returns a tibble with the columns as they appear in all csv files.
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
read_plus <- function(file_list,
                      show_col_types = FALSE) {
  # # testing
  # file_list = file_list
  # filename <- file_list[1]
  # show_col_types = FALSE
  
  df_fin <- tibble()
  for(filename in file_list){
    print(filename)
    curr_df <- read_csv(filename,
                        show_col_types = show_col_types) %>% 
      mutate(filename = gsub("\\.csv", "", basename(filename)))
    df_fin <- rbind(df_fin,curr_df)
  }
  return(df_fin)
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