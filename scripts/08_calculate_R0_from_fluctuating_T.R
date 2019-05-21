# calculate R0 with daily ((nocturnal or diurnal) fluctuating temp

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


in_path <- file.path("output", "termal_response_fits", "informative")

dir_save <- file.path("figures", "trait_R0_relationships")

TS_file_name <- "TS2_DayTemp.rds"


# are you using the cluster? --------------------------------------------------  


if (CLUSTER) {

  obj <- didehpc::queue_didehpc(ctx)
  
} else {
  
  context::context_load(ctx)
  context::parallel_cluster_start(8, ctx)
  
}


# load data -------------------------------------------------------------------


a_samps <- readRDS(file.path(in_path, "a_samps.rds"))
b_samps <- readRDS(file.path(in_path, "b_samps.rds"))
c_samps <- readRDS(file.path(in_path, "c_samps.rds"))
MDR_samps <- readRDS(file.path(in_path, "MDR_samps.rds"))
EFD_samps <- readRDS(file.path(in_path, "EFD_samps.rds")) 
e2a_samps <- readRDS(file.path(in_path, "e2a_samps.rds"))
PDR_samps <- readRDS(file.path(in_path, "PDR_samps.rds"))
lf_samps <- readRDS(file.path(in_path, "lf_DENV_samps.rds"))

foi_data <- read.csv(file.path("output", "extracted_covariates.csv"))

TS2 <- readRDS(file.path("output", "trait_R0_relationships", TS_file_name))


# pre processing --------------------------------------------------------------


no_data <- nrow(foi_data)
n <- dim(a_samps)[1]
thinned <- seq(1, n, by = 5) 
lthin <- length(thinned)

  
# submit all jobs ------------------------------------------------------------- 


if (CLUSTER) {
  
  annual_R0_all_data_points <- queuer::qlapply(
    seq_len(no_data),
    calculate_annual_traits_R0,
    obj,
    lthin,  
    TS2, 
    thinned, 
    a_samps, 
    PDR_samps, 
    MDR_samps, 
    e2a_samps, 
    b_samps, 
    c_samps, 
    lf_samps, 
    EFD_samps)
  
} else {
  
  annual_R0_all_data_points <- loop(seq_len(no_data), 
                                    calculate_annual_traits_R0,
                                    lthin,  
                                    TS2, 
                                    thinned, 
                                    a_samps, 
                                    PDR_samps, 
                                    MDR_samps, 
                                    e2a_samps, 
                                    b_samps, 
                                    c_samps, 
                                    lf_samps, 
                                    EFD_samps,
                                    parallel = TRUE)
  
}
