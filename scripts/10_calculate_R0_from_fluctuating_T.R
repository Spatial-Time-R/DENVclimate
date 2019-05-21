# plot 


# define parameters -----------------------------------------------------------


dir_save <- file.path("figures", "trait_R0_relationships")

covariates <- c("NightTemp_const_term", "DayTemp_const_term")

var <- "pred_R0_1"


# define variables ------------------------------------------------------------


covar <- covariates[1]


# load data -------------------------------------------------------------------


foi_data <- read.csv(file.path("output", "extracted_covariates.csv"))

R0.M <- readRDS(file.path("output", "trait_R0_relationships", paste0(var, "_", covar, "_fluctuating_T.rds")))
  

# pre processing --------------------------------------------------------------


foi_data$pred_R0_1 <- R0.M

x_range <- pretty(foi_data[, "R0_1"])
y_range <- pretty(foi_data[, var])
  
lm <- lm(as.formula(paste0(var, "~ R0_1 - 1")), data = foi_data)
r_sq <- round(summary(lm)$r.squared, 3)


# make plots ------------------------------------------------------------------


dir.create(dir_save, FALSE, TRUE)

png(file.path(dir_save, sprintf("obs_vs_%s_%s%s", var, covar, "_fluctuating_T.png")), 
    width = 9, 
    height = 8, 
    units = "cm",
    pointsize = 12,
    res = 200)

par(mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")

plot(foi_data[, "R0_1"], 
     foi_data[, var], 
     xlim = c(0, max(x_range)),
     ylim = c(0, max(y_range)),
     xlab = "Observations", 
     ylab = "Predictions", 
     pch = 19, 
     cex = 0.5, 
     axes = FALSE)

title(covar, cex.main = 1)
axis(side = 1, at = x_range)
axis(side = 2, at = y_range, las = 2)

abline(reg = lm, col = "red", lwd = 2)
text(10, 6, labels = bquote(R^2 == .(r_sq)), col = "red", lwd = 2)

dev.off()
