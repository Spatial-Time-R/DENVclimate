# plot


# define parameters -----------------------------------------------------------


dir_save <- file.path("figures", "trait_R0_relationships")

covariates <- c("NightTemp_const_term", "DayTemp_const_term")

var <- "pred_R0_1"


# start -----------------------------------------------------------------------


for (i in seq_along(covariates)){

  covar <- covariates[i]


  # load data -------------------------------------------------------------------


  foi_covariates <- readRDS(file.path("output", "extracted_covariates.rds"))

  R0.M <- readRDS(file.path("output", "trait_R0_relationships", paste0(var, "_", covar, "_fluctuating_T.rds")))


  # pre processing --------------------------------------------------------------


  foi_covariates$pred_R0_1 <- R0.M


  # make basic scatter plot and save

  p <- basic_scatter_plot(df = foi_covariates,
                          x = covar,
                          y = var)
  save_plot(p,
            out_pth = dir_save,
            out_fl_nm = sprintf("%s_%s_%s", var, covar, "fluctuating_T.png"),
            wdt = 8,
            hgt = 8)

  # make predictions vs observations plot

  p2 <- preds_vs_obs_scatter_plot(df = foi_covariates,
                                  x = "R0_1",
                                  y = var,
                                  covar = covar)

  save_plot(p2,
            out_pth = dir_save,
            out_fl_nm = sprintf("obs_vs_%s_%s_%s", var, covar, "fluctuating_T.png"),
            wdt = 9,
            hgt = 8)

}
