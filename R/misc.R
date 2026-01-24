#' Create folders
#'
#' Creates folders in folder structure expected by the
#' CompoundVision3D package.
#'
#' @param parent_folder Path to where folders should be created. 
#' If `NULL`, current working directory will be used.
#' Default: `NULL`. 
#'
#' @return Creates the following folders:\cr
#' `./data/`\cr
#' `./data/1_pre_STLs/`\cr
#' `./data/2_STLs/`\cr
#' `./data/3_triangle_centers_and_normals/`\cr
#' `./data/4_local_heights/`\cr
#' `./data/5_rough_clusters/`\cr
#' `./data/6_facet_candidates/`\cr
#' `./data/7_facet_positions/`\cr
#' `./data/8_facets_landmarks_combined/`\cr
#' `./data/9_facets_landmarks_rotated/`\cr
#' `./data/10_facet_infos/`\cr
#' `./data/11_analyses/`
#' @importFrom openxlsx write.xlsx
#' @importFrom tibble tibble
#' @importFrom dplyr add_row
#'
#' @export
#' @examples
#' xxx: add example
#'
set_up_folder_structure <- function(parent_folder = NULL,
                                    verbose = TRUE){
  # # testing
  # parent_folder <- "X:/Pub/2025/Ruehr_Mayfly_eyes/analysis"
  
  if(is.null(parent_folder)){
    parent_folder = getwd()
  }
  
  # Folder definition ---------------------------------------------------------
  # define and create all folders needed in the workflow
  
  # define parent folder for other folders (called data)
  data_folder <- file.path(parent_folder, "data")
  
  # define raw data folder
  raw_data_folder <- "0_raw_data"
  
  # define parent folder for all STLs
  pre_stl_folder <- "1_pre_STLs"
  
  # define folder with STLs to import
  stl_folder <- "2_STLs"
  
  # define folder to save output (triangle centers and normals)
  triangle_centers_and_normals_folder <- "3_triangle_centers_and_normals"
  
  # define folder to store local heights
  local_heights_folder <- "4_local_heights"
  
  # # define folder to store normalized local heights
  # local_heights_normalized_folder <- "4_2_local_heights_normalized"
  
  # define folder to store rough clusters
  rough_clusters_folder <- "5_rough_clusters"
  
  # define folder to store fine clusters
  facet_candidate_folder <- "6_facet_candidates"
  
  # define folder for Blender corrections of fine clusters
  facet_positions_folder <- "7_facet_positions"
  
  # define folder for to save corrected facet positions plus landmarks combined
  facets_landmarks_combined_folder <- "8_facets_landmarks_combined"
  
  # define folder for to save rotated facet positions plus landmarks combined
  facets_landmarks_combined_rotated_folder <- "9_facets_landmarks_rotated"
  
  # define folder for to save infos of rotated facets
  facet_infos_folder <- "10_facet_infos"
  
  # define folder for for further analyses
  analysis_folder <- "11_analyses"
  
  folders_to_create = c(raw_data_folder,
                        pre_stl_folder, stl_folder,
                        triangle_centers_and_normals_folder,
                        local_heights_folder, # local_heights_normalized_folder,
                        rough_clusters_folder, 
                        facet_candidate_folder, facet_positions_folder,
                        facets_landmarks_combined_folder, facets_landmarks_combined_rotated_folder,
                        facet_infos_folder, analysis_folder)
  
  # Folder creation ---------------------------------------------------------
  
  if(dir.exists(data_folder)){
    if(verbose == TRUE) cat("Folder", data_folder, "already exists.\n")
  } else{
    if(verbose == TRUE) cat("Creating folder", data_folder, "...\n")
    dir.create(data_folder)
  }
  # create folders if they do not yet exist
  for(folder in folders_to_create){
    folder_path = file.path(file.path(parent_folder, "data"), folder)
    if(dir.exists(folder_path)){
      if(verbose == TRUE) cat("Folder", folder_path, "already exists.\n")
    } else{
      if(verbose == TRUE) cat("Creating", folder_path, "...\n")
      dir.create(folder_path)
    }
  }
  
  
  raw_data_folders <- list.dirs(file.path(file.path(parent_folder, "data"), raw_data_folder),
                                full.names = TRUE)[-1]
  if(length(raw_data_folders) > 0){
    CV_numbers <- c()
    if (verbose == TRUE) cat("Creating .\\1_pre_STLs\\3D_Slicer for each folder in .\\0_raw_data...\n")
    for(folder in raw_data_folders){
      dataset_preSTL_folder = file.path(gsub(raw_data_folder, pre_stl_folder, folder))
      slicer_folder = file.path(dataset_preSTL_folder, "3D_Slicer")
      if(dir.exists(dataset_preSTL_folder)){
        if(verbose == TRUE) cat("Folder", dataset_preSTL_folder, "already exists.\n")
      } else{
        if(verbose == TRUE) cat("Creating", dataset_preSTL_folder, "...\n")
        dir.create(dataset_preSTL_folder)
      }
      if(dir.exists(slicer_folder)){
        if(verbose == TRUE) cat("Folder", slicer_folder, "already exists.\n")
      } else{
        if(verbose == TRUE) cat("Creating", slicer_folder, "...\n")
        dir.create(slicer_folder)
      }
      
      # extract CV numbers
      curr_CV_number = gsub("^CV(\\d{4})_.+", "\\1", basename(dataset_preSTL_folder))
      CV_numbers <- c(CV_numbers, curr_CV_number)
    }
  }
  
  # excel sheets
  excel_file_path <- file.path(data_folder, "CompoundVision3D_EyeNotes.xlsx")
  if(file.exists(excel_file_path)){
    if(verbose == TRUE) cat(excel_file_path, "already exists.\n")
  } else {
    if(verbose == TRUE) cat("Creating", excel_file_path, "...\n")
    
    # CompoundVision3D_EyeNotes
    # library(openxlsx)
    # library(tibble)
    # library(dplyr)
    CompoundVision3D_EyeNotes <- tibble(CV = character(),
                                        class = character(),
                                        order = character(),
                                        family = character(),
                                        genus = character(),
                                        species = character(),
                                        eye = character(),
                                        side = character(),
                                        facet_size_estimate = character())
    print(CompoundVision3D_EyeNotes)
    print(CV_numbers)
    print(length(CV_numbers))
    if(length(CV_numbers > 0)){
      for(curr_CV in CV_numbers){
        CompoundVision3D_EyeNotes <- CompoundVision3D_EyeNotes %>% 
          add_row(CV = curr_CV,
                  eye = "1") %>% 
          add_row(CV = curr_CV,
                  eye = "2")
      }
    }
    print(CompoundVision3D_EyeNotes)
    
    write.xlsx(CompoundVision3D_EyeNotes,
               rowNames = FALSE,
               showNA = FALSE,
               sheetName = "EyeNotes",
               file =  excel_file_path)
  }
}
