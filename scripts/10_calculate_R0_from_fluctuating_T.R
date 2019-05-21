# load results and save 

options(didehpc.cluster = "fi--didemrchnb")

CLUSTER <- TRUE

my_resources <- c(
  file.path("R", "Mordecai_et_al_2017_code", "temp_functions_all.R"),
  file.path("R", "R0_function.R"),
  file.path("R", "calculate_annual_traits_R0_FT.R"),
  file.path("R", "utility_functions.R"))

context::context_log_start()
ctx <- context::context_save(path = "context",
                             sources = my_resources)


# define parameters -----------------------------------------------------------


dir_save <- file.path("output", "trait_R0_relationships")

covariates <- c("NightTemp_const_term", "DayTemp_const_term")

var <- "pred_R0_1"


# define variables ------------------------------------------------------------


covar <- covariates[1]


# load data -------------------------------------------------------------------


foi_data <- read.csv(file.path("output", "extracted_covariates.csv"))


# are you using the cluster? --------------------------------------------------  


if (CLUSTER) {
  
  obj <- didehpc::queue_didehpc(ctx)
  
} else {
  
  context::context_load(ctx)
  context::parallel_cluster_start(8, ctx)
  
}


# get the results -------------------------------------------------------------


all_bundles <- obj$task_bundle_info()

id <- all_bundles[nrow(all_bundles), "name"]

task_obj <- obj$task_bundle_get(id)

all_results <- task_obj$results()

all_results_mat <- do.call("rbind", all_results)

R0.M <- rowMeans(all_results_mat)


# save ------------------------------------------------------------------------


write_out_rds(R0.M, dir_save, paste0(var, "_", covar, "_fluctuating_T.rds"))
