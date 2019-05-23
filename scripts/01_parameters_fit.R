### This file will process the temperature data and create the individual parameter
### samples for all the traits that are included in the DENV R0 model.


# define parameters -----------------------------------------------------------


n.chains <- 5
n.adapt <- 5000
n.samps <- 5000
Temps <- seq(0, 50, by = 0.1)
priors1 <- list()
priors1$names <- c("T0", "Tm", "c","tau")
priors1$fun <- c("uniform", "uniform", "gamma","gamma")
priors1$hyper <- matrix(NA, ncol = 4, nrow = 3)
priors1$hyper[, 1] <- c(0, 24, NA)
priors1$hyper[, 2] <- c(25, 45, NA)
priors1$hyper[, 3] <- c(1, 10, NA)
priors1$hyper[, 4] <- c(0.0001, 0.0001, NA)
thinned_vec <- seq(1, n.samps, length = 1000)

jags_briere_path <- system.file("extdata", "jags-briere.bug", package = "DENVclimate")
jags_briere_EFD_path <- system.file("extdata", "jags-briere-EFD.bug", package = "DENVclimate")
jags_briere_trunc_b_path <- system.file("extdata", "jags-briere-trunc-b.bug", package = "DENVclimate")
jags_briere_trunc_path <- system.file("extdata", "jags-briere-trunc.bug", package = "DENVclimate")
jags_quad_neg_path <- system.file("extdata", "jags-quad-neg.bug", package = "DENVclimate")

out_path <- file.path("output", "termal_response_fits", "uninformative")


# load data -------------------------------------------------------------------


# Exclude the Focks & Barrera 2006 data because they're from a model
data.all <- subset(trait_data,
                   ref != as.character("Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"))



# -----------------------------------------------------------------------------
#
# GCR
#
# -----------------------------------------------------------------------------



data <- data.all[which(data.all$trait.name == "GCR"), ]

# plot(trait ~ T, data = data)

jags <- rjags::jags.model(jags_briere_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 31, T0 = 5, c = 0.00007),
                          n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- rjags::coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnostic information.

# plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
a.samps <- samps

# plot priors with histogram of samples
# plot.hists(a.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

# calculate the response
out1 <- make.sims.temp.resp(sim = "briere", a.samps, Temps, thinned = thinned_vec)

# calculate the 95% inner quantile/HPD
q1 <- temp.sim.quants(out1$fits, length(Temps))

# plot the data with the fits/quantiles
# plot(data$T, data$trait,
#      xlim = c(10, 40),
#      ylim = c(0, 0.42),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out1$fits, q = q1)



# -----------------------------------------------------------------------------
#
# EFD
#
# -----------------------------------------------------------------------------



# number of eggs laid per female per day

data <- data.all[which(data.all$trait.name == "EFD"), ]

# plot(trait ~ T, data = data)

nlfit <- nls(trait ~ briere(T, c, Tm, T0),
             data = data,
             start = list(c = 0.1, Tm = 35, T0 = 10),
             algorithm = "port",
             lower = c(0, max(data$T), 0))

# lines(Temps, briere(Temps, coef(nlfit)[1], coef(nlfit)[2], coef(nlfit)[3]))

jags <- rjags::jags.model(jags_briere_EFD_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 34, T0 = 15, c = 0.007),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

#plot(coda.samps)

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
EFD.samps <- samps

# plot.hists(EFD.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out2 <- make.sims.temp.resp(sim = "briere", EFD.samps, Temps, thinned = thinned_vec)

q2 <- temp.sim.quants(out2$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out2$fits, q = q2)



# -----------------------------------------------------------------------------
#
# b
#
# -----------------------------------------------------------------------------



## The next trait we'll be choosing is the vector competence, b*c.
## Due to the data that was collected it will be necessary
## to decompose it into its two parts, b and c.

# Choose b, the probability a human will be bitten and infected by
# an infectious mosquito (ie. transmission).

data <- data.all[which(data.all$trait.name == "b"),]

# remove Lambrechts data since they are from other mosquitoes/flaviviruses
data <- subset(data, ref!="Lambrects_et_al_2011_PNAS")

# plot(trait ~ T, data = data, xlim = c(15, 40), ylim = c(0, 1))
# points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
# points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)

jags <- rjags::jags.model(jags_briere_trunc_b_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 31, T0 = 5, c = 0.00007),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# plot(coda.samps)

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
b.samps <- samps

# plot.hists(b.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out3 <- make.sims.temp.resp(sim = "briere.trunc", b.samps, Temps, thinned = thinned_vec)

q3 <- temp.sim.quants(out3$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      ylim = c(0, 1),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out3$fits, q = q3)
# points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
# points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)



# -----------------------------------------------------------------------------
#
# c
#
# -----------------------------------------------------------------------------



# probability that a mosquito becomes infected after biting an infectious human (ie. infection).

data <- data.all[which(data.all$trait.name == "c"), ]

# remove Lambrechts data since they are from other mosquitoes/flaviviruses
data <- subset(data, ref!="Lambrects_et_al_2011_PNAS")
data <- subset(data, ref!="Alto&Bettinardi_2013_AJTMH")

# plot(trait ~ T, data = data)
# points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
# points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)

jags <- rjags::jags.model(jags_briere_trunc_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 31, T0 = 5, c = 0.00007),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# plot(coda.samps)

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
c.samps <- samps

# plot.hists(c.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out4 <- make.sims.temp.resp(sim = "briere.trunc", c.samps, Temps, thinned = thinned_vec)

q4 <- temp.sim.quants(out4$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      ylim = c(0, 1),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out4$fits, q = q4)

### Older version of the code, DENV_IndividualParameterFitCode.R
### fits a quadratic function to c to compare. Omitted here.



# -----------------------------------------------------------------------------
#
# MDR
#
# -----------------------------------------------------------------------------



# mean development time for a mosquito.

data <- data.all[which(data.all$trait.name == "MDR"), ]

# plot(trait ~ T, data = data)

jags <- rjags::jags.model(jags_briere_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 31, T0 = 5, c = 0.00007),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# plot(coda.samps)

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
MDR.samps <- samps

# plot.hists(MDR.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out5 <- make.sims.temp.resp(sim = "briere", MDR.samps, Temps, thinned = thinned_vec)

q5 <- temp.sim.quants(out5$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out5$fits, q = q5)



# -----------------------------------------------------------------------------
#
# pEA
#
# -----------------------------------------------------------------------------



# probability a mosquito will survive from hatching to maturation.

data <- data.all[which(data.all$trait.name == "pEA"), ]

# plot(trait ~ T, data = data)

jags <- rjags::jags.model(jags_quad_neg_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                          n.chains = n.chains,
                          inits = list(T0=5, Tm=33, n.qd=0.005),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('T0','Tm', 'qd'), n.samps)

# plot(coda.samps)

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, samp.lims = c(1, n.samps), sig = FALSE)
samps.q$n.qd <- samps.q$qd
e2a.samps <- samps.q

# plot.hists(e2a.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out6 <- make.sims.temp.resp(sim = "quad.trunc", e2a.samps, Temps, thinned = thinned_vec)

q6 <- temp.sim.quants(out6$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out6$fits, q = q6)



# -----------------------------------------------------------------------------
#
# p
#
# -----------------------------------------------------------------------------



# daily probability that an adult mosquito survives.

data.days <- data.all[which(data.all$trait.name == "p/days"), ]
data.1mu <- data.all[which(data.all$trait.name == "1/mu"), ]

data.days$trait <- 1/data.days$trait

data <- rbind(data.days, data.1mu)

# plot(trait ~ T, data = data)
# points(trait ~ T, data = subset(data, trait.name=="p/days"), col=2, pch=16)
# points(trait ~ T, data = subset(data, trait.name=="1/mu"), col=4, pch=16)

jags <- rjags::jags.model(jags_quad_neg_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                          n.chains = n.chains,
                          inits = list(T0=5, Tm=33, n.qd=0.005),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('T0','Tm', 'qd', 'sigma'), n.samps)

# plot(coda.samps, ask = TRUE)

samps.q <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE)

samps.q$n.qd <- samps.q$qd
samps.q$tau <- 1/samps.q$sigma
lf.DENV.samps <- samps.q

priors1 <- list()
priors1$names <- c("inter", "n.slope", "qd","tau")
priors1$fun <- c("uniform", "uniform", "gamma","gamma")
priors1$hyper <- matrix(NA, ncol = 4, nrow = 3)
priors1$hyper[, 1] <- c(0, 24, NA)
priors1$hyper[, 2] <- c(25, 45, NA)
priors1$hyper[, 3] <- c(1, 1, NA)
priors1$hyper[, 4] <- c(0.0001, 0.0001, NA)

# plot.hists(cbind(lf.DENV.samps[,c(1,2)], -lf.DENV.samps[,3], lf.DENV.samps[,4]), my.par=c(2,2), n.hists=4, priors=priors1)

out7 <- make.sims.temp.resp(sim = "quad", lf.DENV.samps, Temps, thinned = seq(1, 25000, by = 5), trunc.num = 0.0001)

q7 <- temp.sim.quants(out7$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 40),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = "Lifespan (days)",
#      cex = 2)
# add.sim.lines(Temps, sim.data = out7$fits, q = q7)



# -----------------------------------------------------------------------------
#
# PDR
#
# -----------------------------------------------------------------------------



# parasite/virus development rate.

data1 <- data.all[which(data.all$trait.name=="PDR"),]
data2 <- data.all[which(data.all$trait.name=="EIP"),]
data2$trait <- 1/data2$trait
data <- rbind(data1, data2)

# Plot the data to see which function, Briere or Quadratic, is a more suitable fit for
# the data.

# plot(trait ~ T, data = data)
# points(trait ~ T, data = subset(data, ref=="Davis_1932_AmJEpidemiology"), col=2, pch=16)
# points(trait ~ T, data = subset(data, ref=="Focks_et_al_1995_AJTMH"), col=3, pch=16)
# points(trait ~ T, data = subset(data, ref=="McLean_et_al_1974_CanJMicobiol"), col=4, pch=16)
# points(trait ~ T, data = subset(data, ref=="McLeah_et_al_1975_MosquitoNews"), col=5, pch=16)
# points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=6, pch=16)
# legend('topleft', legend=c("Davis", "Focks", "McLean 74", "McLean 75", "Watts", "Carrington"), col=c(2:6, 1), pch=c(rep(16, 5), 1))

jags <- rjags::jags.model(jags_briere_path,
                          data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                          n.chains = n.chains,
                          inits = list(Tm = 38, T0 = 5, c = 0.00007),
                          n.adapt = n.adapt)

coda.samps <- rjags::coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# plot(coda.samps)

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
PDR.samps <-  samps

# plot.hists(PDR.samps[, c(1:3, 4)], my.par = c(2, 2), n.hists = 4, priors = priors1)

out8 <- make.sims.temp.resp(sim = "briere", PDR.samps, Temps, thinned = thinned_vec)

q8 <- temp.sim.quants(out8$fits, length(Temps))

# plot(data$T,
#      data$trait,
#      xlim = c(10, 45),
#      pch = 20,
#      xlab = "Temperature (C)",
#      ylab = data$trait.name[1],
#      cex = 2)
# add.sim.lines(Temps, sim.data = out7$fits, q = q7)
# points(trait ~ T, data = subset(data, ref=="Davis_1932_AmJEpidemiology"), col=2, pch=16)
# points(trait ~ T, data = subset(data, ref=="Focks_et_al_1995_AJTMH"), col=3, pch=16)
# points(trait ~ T, data = subset(data, ref=="McLean_et_al_1974_CanJMicobiol"), col=4, pch=16)
# points(trait ~ T, data = subset(data, ref=="McLeah_et_al_1975_MosquitoNews"), col=5, pch=16)
# points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=6, pch=16)
# legend('topleft', legend=c("Davis", "Focks", "McLean 74", "McLean 75", "Watts", "Carrington"), col=c(2:6, 1), pch=c(rep(16, 5), 1))


# save ------------------------------------------------------------------------


# save(a.samps, b.samps, c.samps, MDR.samps, EFD.samps, e2a.samps, PDR.samps, lf.DENV.samps,
#      file = file.path("output", "Aegypti_DENV_ParameterFits_2016-03-30.Rsave"))
#
# save(PDR.samps, b.samps, c.samps,
#      file = file.path("output", "Aegypti_DENV_b_c_PDR_samps_2016-03-30.Rsave"))

write_out_rds(a.samps, out_path, "a_samps.rds")
write_out_rds(b.samps, out_path, "b_samps.rds")
write_out_rds(c.samps, out_path, "c_samps.rds")
write_out_rds(MDR.samps, out_path, "MDR_samps.rds")
write_out_rds(EFD.samps, out_path, "EFD_samps.rds")
write_out_rds(e2a.samps, out_path, "e2a_samps.rds")
write_out_rds(PDR.samps, out_path, "PDR_samps.rds")
write_out_rds(PDR.samps, out_path, "lf.DENV.samps")


# =============================================================================
#                  plot data and fits of all traits together
# =============================================================================


selected_traits <- c("GCR", "EFD", "b", "c", "MDR", "pEA", "p", "EIP")

# subset original dataset
data.all <- subset(data.all, !(trait.name == "b" & ref == "Lambrects_et_al_2011_PNAS"))
data.all <- subset(data.all, !(trait.name == "c" & (ref == "Lambrects_et_al_2011_PNAS" | ref == "Alto&Bettinardi_2013_AJTMH")))

my_data <- data.all[data.all$trait.name %in% selected_traits[selected_traits != "EIP" & selected_traits != "p"], ]

data.days <- data.all[which(data.all$trait.name == "p/days"), ]
data.1mu <- data.all[which(data.all$trait.name == "1/mu"), ]
data.days$trait <- 1/data.days$trait
sub_data_1 <- rbind(data.days, data.1mu)
sub_data_1[which(sub_data_1[, "trait.name"] == "p/days"), "trait.name"] <- "p"
sub_data_1[which(sub_data_1[, "trait.name"] == "1/mu"), "trait.name"] <- "p"

data1 <- data.all[which(data.all$trait.name=="PDR"),]
data2 <- data.all[which(data.all$trait.name=="EIP"),]
data2$trait <- 1/data2$trait
sub_data_2 <- rbind(data1, data2)
sub_data_2[which(sub_data_2[, "trait.name"] == "PDR"), "trait.name"] <- "EIP"

my_data <- rbind(my_data, sub_data_1, sub_data_2)

my_data$trait.name <- droplevels(my_data$trait.name)

my_data$trait.name <- factor(my_data$trait.name, levels = selected_traits)

all_T <- rep(out1$T, length(selected_traits))
all_out <- list(out1$fits, out2$fits, out3$fits, out4$fits, out5$fits, out6$fits, out7$fits, out8$fits)
all_q <- list(q1, q2, q3, q4, q5, q6, q7, q8)
all_q_t <- lapply(all_q, t)
all_out_m <- lapply(all_out, rowMeans)
all_q_t_mat <- do.call("rbind", all_q_t)
colnames(all_q_t_mat) <- c("q1", "q2")
trait_factor <- rep(selected_traits, each = nrow(all_q_t[[1]]))

all_fits <- cbind(temp = all_T,
                  mean = unlist(all_out_m),
                  as.data.frame(all_q_t_mat),
                  trait.name = trait_factor)


# order factor ----------------------------------------------------------------


ordered_selected_traits <- c("GCR", "EFD", "pEA", "MDR", "p", "b", "c", "EIP")

my_data$trait.name <- factor(my_data$trait.name, levels = ordered_selected_traits)

all_fits$trait.name <- factor(all_fits$trait.name, levels = ordered_selected_traits)


# plot and save ---------------------------------------------------------------


plot_thermal_responses_all_data(my_data,
                                all_fits,
                                "figures",
                                "all_thermal_responses_uninf_priors")


# get mean and se -------------------------------------------------------------

# calculate mean and se of the data for each T value

my_data$T <- round(my_data$T)

my_data_av <- mean_trait(my_data)

my_data_av$se <- my_data_av$sd / sqrt(my_data_av$n)
my_data_av[is.na(my_data_av$se), "se"] <- 0


# plot and save ---------------------------------------------------------------

plot_thermal_responses_mean_data(my_data_av,
                                 all_fits,
                                 "figures",
                                 "all_thermal_responses_uninf_priors_2")
