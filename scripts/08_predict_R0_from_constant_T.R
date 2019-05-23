
# define parameters -----------------------------------------------------------


in_path <- file.path("output", "termal_response_fits", "informative")

covariates <- c("DayTemp_const_term", "NightTemp_const_term")

response <- "pred_R0_1"

dir_save <- file.path("figures", "trait_R0_relationships")


# load data -------------------------------------------------------------------


a_samps <- readRDS(file.path(in_path, "a_samps.rds"))
b_samps <- readRDS(file.path(in_path, "b_samps.rds"))
c_samps <- readRDS(file.path(in_path, "c_samps.rds"))
MDR_samps <- readRDS(file.path(in_path, "MDR_samps.rds"))
EFD_samps <- readRDS(file.path(in_path, "EFD_samps.rds"))
e2a_samps <- readRDS(file.path(in_path, "e2a_samps.rds"))
PDR_samps <- readRDS(file.path(in_path, "PDR_samps.rds"))
lf_samps <- readRDS(file.path(in_path, "lf_DENV_samps.rds"))

foi_covariates <- readRDS(file.path("output", "extracted_covariates.rds"))


# calculate R0 with constant temp ---------------------------------------------


for (i in seq_along(covariates)){

  covar <- covariates[i]

  temp <- foi_covariates[, covar]

  t <- length(temp)

  # Length of samples, assuming they're all the same length
  n <- dim(a_samps)[1]
  thinned <- seq(1, n, by = 5)
  lthin <- length(thinned)

  ### Calculate R0 and each trait across the thinned posterior samples
  R0 <- matrix(NA, t, lthin)
  a <- b <- c <- PDR <- MDR <- EFD <- e2a <- lf <- matrix(NA, t, lthin)

  for (j in seq_len(lthin)) {

    # if(j %% 50 == 0) cat("iteration =", j, "\n")

    # calculate parameter trajectories
    i <- thinned[j]
    a[, j] <- briere(temp, a_samps[i,3], a_samps[i,2], a_samps[i,1])
    PDR[, j] <- briere(temp, PDR_samps[i,3], PDR_samps[i,2], PDR_samps[i,1])
    MDR[, j] <- briere(temp, MDR_samps[i,3], MDR_samps[i,2], MDR_samps[i,1])
    EFD[, j] <- briere(temp, EFD_samps[i,3], EFD_samps[i,2], EFD_samps[i,1])
    e2a[, j] <- quad.2.trunc(temp, e2a_samps[i,1], e2a_samps[i,2], e2a_samps[i,3])
    b[, j] <- briere.trunc(temp, b_samps[i,3], b_samps[i,2], b_samps[i,1])
    c[, j] <- briere.trunc(temp, c_samps[i,3], c_samps[i,2], c_samps[i,1])
    lf[, j] <- quad.2(temp, lf_samps[i,1], lf_samps[i,2], lf_samps[i,3])

    # Calculate R0 equation
    R0[, j] <- myR0(a[, j], b[, j], c[, j], PDR[, j], MDR[, j], EFD[, j], e2a[, j], lf[, j])

  }

  R0.M <- rowMeans(R0)
  foi_covariates$pred_R0_1 <- R0.M

  # make plot and save

  p <- basic_scatter_plot(df = foi_covariates,
                          x = covar,
                          y = response)

  save_plot(p,
            out_pth = dir_save,
            out_fl_nm = paste0(response, "_", covar),
            wdt = 8,
            hgt = 8)

}
