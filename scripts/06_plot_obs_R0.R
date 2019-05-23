
# plot R0 vs mean annual day and night temperatures


# define parameters -----------------------------------------------------------


response <- "R0_1"

covariates <- c("DayTemp_const_term", "NightTemp_const_term")

dir_save <- file.path("figures", "trait_R0_relationships")


# load data -------------------------------------------------------------------


foi_covariates <- readRDS(file.path("output", "extracted_covariates.rds"))


# make plots ------------------------------------------------------------------


for (i in seq_along(covariates)){

  dir.create(dir_save, FALSE, TRUE)

  covar <- covariates[i]

  png(file.path(dir_save, paste0(response, "_", covar,".png")),
      width = 8,
      height = 8,
      units = "cm",
      pointsize = 12,
      res = 200)

  par(mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))

  plot(foi_covariates[, covar],
       foi_covariates[, response],
       xlab = covar,
       ylab = response,
       pch = 19,
       cex = 0.5)

  j <- order(foi_covariates[, covar])
  l_1 <- loess(as.formula(paste0(response, "~", covar)), data = foi_covariates)
  lines(foi_covariates[, covar][j], l_1$fitted[j], col = "red", lwd = 3)

  dev.off()

}
