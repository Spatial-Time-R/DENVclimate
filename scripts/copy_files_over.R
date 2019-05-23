# copy over files from shared drive

root <- file.path("Y:", "DENV_risk_maps")

paths_from <- list(
  file.path("data", 
            "foi", 
            "FOI_estimates_lon_lat.csv"),
  file.path("output", 
            "env_variables", 
            "All_adm1_env_var.csv"),
  file.path("output", 
            "datasets", 
            "country_age_structure.csv"),
  file.path("output",
            "datasets",
            "country_age_structure_mean.csv"),
  file.path("output", 
            "env_variables", 
            "All_adm1_env_var.csv"),
  file.path("output", 
            "datasets", 
            "pseudo_absence_points_2.csv"),
  file.path("data", 
            "env_variables", 
            "plus60minus60_tiles.csv"),
  file.path("output", 
            "datasets", 
            "NA_pixel_tiles_20km.txt"),
  file.path("output", 
            "variable_selection",
            "stepwise_v3",
            "predictor_rank.csv"),
  file.path("output", 
            "env_variables", 
            "all_squares_env_var_0_1667_deg.rds"),
  file.path("output", 
            "shapefiles",
            "gadm28_adm1_eras.shp"),
  file.path("output", 
            "shapefiles",
            "gadm28_adm1_dengue.shp"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_I_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_C_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_HC_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_R0_1_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_R0_2_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_R0_3_lookup_tables.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_C_lookup_tables_fixed_params.rds"),
  file.path("output", 
            "predictions_world", 
            "FOI_to_HC_lookup_tables_fixed_params.rds"))

for (i in seq_along(paths_from)){
  
  path_i <- paths_from[[i]]
  
  message(path_i)
  
  path_from_i <- file.path(root, path_i) 
  
  # match everything before the last occurrence of /
  path_to_i <- sub("/([^/]*)$", "", path_i)
  
  if(!dir.exists(path_to_i)) dir.create(path_to_i, FALSE, TRUE)
  
  file.copy(path_from_i, path_to_i, overwrite = FALSE)

  last_three_digits <- sub("^([^.]*).", "", path_i)
  
  if(last_three_digits == "shp"){
    
    evr_but_last_three_digits <- sub(".([^.]*)$", "", path_i)
    
    path_from_i_1 <- file.path(root, paste0(evr_but_last_three_digits, ".dbf"))  
    file.copy(path_from_i_1, path_to_i, overwrite = FALSE)
    
    path_from_i_2 <- file.path(root, paste0(evr_but_last_three_digits, ".prj"))  
    file.copy(path_from_i_2, path_to_i, overwrite = FALSE)
    
    path_from_i_3 <- file.path(root, paste0(evr_but_last_three_digits, ".shx"))  
    file.copy(path_from_i_3, path_to_i, overwrite = FALSE)
    
  }
  
}
