# create data

trait_data <- read.csv(file.path("data-raw",
                                 "aegyptiDENVmodelTempData_2016-03-30.csv"),
                       header = TRUE)

aedes_aegypti_priors <- read.csv(file.path("data-raw", "Aedes_prior_data.csv"),
                                 header = TRUE)

aa_EIP_priors <- read.csv(file.path("data-raw", "EIP_priors_2015-12-04.csv"),
                          header = TRUE)

foi_data <- read.csv(file.path("data-raw",
                               "All_FOI_estimates_and_predictors.csv"),
                     stringsAsFactors = FALSE)

predictor_rank <- read.csv(file.path("data-raw",
                                     "predictor_rank.csv"),
                           stringsAsFactors = FALSE)

usethis::use_data(trait_data,
                  aedes_aegypti_priors,
                  aa_EIP_priors,
                  foi_data,
                  predictor_rank,
                  overwrite = TRUE)
