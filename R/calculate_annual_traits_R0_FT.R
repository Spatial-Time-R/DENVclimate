calculate_annual_traits_R0 <- function(i,
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
                                       EFD_samps) {

  cat("data point =", i, "\n")

  # Grab the vector of daily temps
  dtemps <- TS2[i, ]

  out <- c()

  # Loop through each iteration in the thinned MCMC chain
  for(j in seq_len(lthin)) {

    # Calculate daily trait and R0 values at different daily temperatures,
    # for each set of function parameters (5000 parameter sets)

    k <- thinned[j]

    if(j %% 500 == 0) cat("iteration =", j, "\n")

    a_daily <- briere(dtemps, a_samps[k, 3], a_samps[k, 2], a_samps[k, 1])
    PDR_daily <- briere(dtemps, PDR_samps[k, 3], PDR_samps[k, 2], PDR_samps[k, 1])
    MDR_daily <- briere(dtemps, MDR_samps[k, 3], MDR_samps[k, 2], MDR_samps[k, 1])
    e2a_daily <- quad.2.trunc(dtemps, e2a_samps[k, 1], e2a_samps[k, 2], e2a_samps[k, 3])
    b_daily <- briere.trunc(dtemps, b_samps[k, 3], b_samps[k, 2], b_samps[k, 1])
    c_daily <- briere.trunc(dtemps, c_samps[k, 3], c_samps[k, 2], c_samps[k, 1])
    lf_daily <- quad.2(dtemps, lf_samps[k, 1], lf_samps[k, 2], lf_samps[k, 3])
    EFD_daily <- briere(dtemps, EFD_samps[k, 3], EFD_samps[k, 2], EFD_samps[k, 1])

    R0_FT_daily <- myR0(a_daily,
                        b_daily,
                        c_daily,
                        PDR_daily,
                        MDR_daily,
                        EFD_daily,
                        e2a_daily,
                        lf_daily)

    # calculate average of daily R0 across the year

    R0_FT_mean <- mean(R0_FT_daily)

    out[j] <- R0_FT_mean

  }

  out

}
