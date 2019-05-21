# create data

trait_data <- read.csv("aegyptiDENVmodelTempData_2016-03-30.csv",
                       header = TRUE)

aedes_aegypti_priors <- read.csv("Aedes_prior_data.csv",
                                 header = TRUE)

aa_EIP_priors <- read.csv("EIP_priors_2015-12-04.csv",
                          header = TRUE)

usethis::use_data(trait_data,
                  aedes_aegypti_priors,
                  aa_EIP_priors,
                  overwrite = TRUE)
