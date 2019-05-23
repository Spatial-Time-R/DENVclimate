# copy over files from shared drive

root <- file.path("Y:", "DENVclimate")

paths_from <- list(
  file.path("data",
            "predictors",
            "env_vars_20km.rds"))

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
