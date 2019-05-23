# Reconstructing time serie from Fourier Transform terms


# define parameters -----------------------------------------------------------


in_path <- file.path("output", "termal_response_fits")

fig_dir_save <- file.path("figures", "trait_R0_relationships")

tab_dir_save <- file.path("output", "trait_R0_relationships")

timepoints <- seq_len(365) / 365

FT_var_names <- c("NightTemp", "DayTemp")

y_lab_tags <- c("Nocturnal", "Diurnal")


# load data -------------------------------------------------------------------


foi_covariates <- readRDS(file.path("output", "extracted_covariates.rds"))


# define variables ------------------------------------------------------------


no_data <- nrow(foi_covariates)


# reconstructing time serie from FT terms -------------------------------------


for (i in seq_along(FT_var_names)){

  FT_var_name <- FT_var_names[i]

  y_lab_tag <- y_lab_tags[i]

  const_term_var_name <- paste0(FT_var_name, "_const_term")
  Re0_var_name <- paste0(FT_var_name, "_Re0")
  Im0_var_name <- paste0(FT_var_name, "_Im0")
  Re1_var_name <- paste0(FT_var_name, "_Re1")
  Im1_var_name <- paste0(FT_var_name, "_Im1")

  TS0 <- matrix(foi_covariates[, const_term_var_name], nrow = no_data, ncol = length(timepoints))

  TS1 <- TS0 + (
    outer(foi_covariates[, Re0_var_name], cos(2 * pi * timepoints * 1)) +
      outer(foi_covariates[, Im0_var_name], sin(2 * pi * timepoints * 1))
  )

  TS2 <- TS1 + (
    outer(foi_covariates[, Re1_var_name], cos(2 * pi * timepoints * 2)) +
      outer(foi_covariates[, Im1_var_name], sin(2 * pi * timepoints * 2))
  )

  write_out_rds(TS2, tab_dir_save, paste0("TS2_", FT_var_name, ".rds"))


  # plotting at one selected location -------------------------------------------


  j <- 1

  y_range <- pretty(c(TS0[j, ], TS1[j, ], TS2[j, ]))

  png(file.path(fig_dir_save, paste0("reconstructed_annual_TS_at_one_location_", FT_var_name, ".png")),
      width = 12,
      height = 10,
      units = "cm",
      pointsize = 12,
      res = 200)

  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

  plot(timepoints, TS0[j,],
       type = "l",
       lwd = 2,
       xlim = c(0, 1),
       ylim = c(min(y_range), max(y_range)),
       xlab = "time throughout the year",
       ylab = bquote(.(y_lab_tag) ~ "Temperature (" ~ degree ~ C ~ ")"),
       main = "",
       axes = FALSE)

  lines(timepoints, TS1[j, ], col = 2, lwd = 2)
  lines(timepoints, TS2[j, ], col = 3, lwd = 3)
  axis(side = 1, at = seq(0, 1, by = 0.2))
  axis(side = 2, at = y_range, las = 2)

  legend("bottomleft",
         col = 1:3,
         lwd = 2,
         legend = c("TS0", "TS1", "TS2"))

  title(main = paste0("Reconstructed time series \nfor longitude ",
                      foi_covariates$longitude[j],
                      " and latitude ",
                      foi_covariates$latitude[j]))

  dev.off()

}
