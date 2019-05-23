# Extract covariates values at data point locations
# Use the 20 km square covariate values for serology data
# and the admin unit pop weighted average for case report data


# define parameters -----------------------------------------------------------


parameters <- list(
  resample_grid_size = 20,
  no_predictors = 26)

grp_flds <- c("ID_0", "ID_1", "data_id")

res <- (1 / 120) * parameters$resample_grid_size

year.i <- 2007

year.f <- 2014

ppyear <- 64

base_info <- c("data_id", "type", "date", "longitude", "latitude", "country", "ISO", "ID_0", "ID_1", "FOI", "R0_1", "R0_2", "R0_3")

out_path <- "output"

out_name <- "extracted_covariates.rds"


# load data -------------------------------------------------------------------


pxl_data <- readRDS(file.path("data",
                              "predictors",
                              "env_vars_20km.rds"))


# pre processing --------------------------------------------------------------


number_of_predictors <- parameters$no_predictors

my_predictors <- predictor_rank$name[1:number_of_predictors]

pxl_data <- inner_join(pxl_data, foi_data[, c(grp_flds, "type")])


# fix serology new_weights ----------------------------------------------------


pxl_data[pxl_data$type == "serology", "new_weight"] <- 0

sero_points <- foi_data[foi_data$type == "serology", ]

pxl_data$lat.int <- round(pxl_data$latitude / res)
pxl_data$long.int <- round(pxl_data$longitude / res)

sero_points$lat.int <- round(sero_points$latitude / res)
sero_points$long.int <- round(sero_points$longitude / res)

sero_points$cell <- 0
sero_points$no_square <- 0

for (i in seq_len(nrow(sero_points))){

  sero_long <- sero_points[i, "long.int"]
  sero_lat <- sero_points[i, "lat.int"]

  matches <- pxl_data$type == "serology" & pxl_data$lat.int == sero_lat & pxl_data$long.int == sero_long

  if(sum(matches) != 0){

    message(i)

    cell_id <- which(matches == TRUE)[1]
    sero_points[i, "cell"] <- cell_id
    pxl_data[cell_id, "new_weight"] <- 1

  } else {

    sero_points[i, "no_square"] <- 1

  }

}

missing_square <- sero_points[sero_points$no_square == 1, ]

sero_covariates <- pxl_data[sero_points$cell, my_predictors]

caseReport_covariates <- foi_data[foi_data$type == "caseReport", my_predictors]

sero_data <- cbind(sero_points[,base_info], sero_covariates)

caseReport_data <- cbind(foi_data[foi_data$type == "caseReport", base_info], caseReport_covariates)

all_data_covariates <- rbind(sero_data, caseReport_data)


# scale covariate values before using R0 model --------------------------------


for (i in seq_along(my_predictors)){

  my_pred <- my_predictors[i]

  scale <- get_covariate_scaling_factor(my_pred, ppyear, year.f, year.i)

  message(scale)

  all_data_covariates[, my_pred] <- all_data_covariates[, my_pred] / scale

}


# save ------------------------------------------------------------------------


write_out_rds(all_data_covariates, out_path, out_name)
