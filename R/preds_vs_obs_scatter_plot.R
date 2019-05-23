

#------------------------------------------------------------------------------
#' preds_vs_obs_scatter_plot
#'
#' \code{preds_vs_obs_scatter_plot} makes a simple scatter plot of predicted vs
#'  oserved data points.
#'
#' @param df dataframe with the data to plot.
#' @param x variable to plot on the x axis.
#' @param y variable to plot on the y axis.
#' @param covar type of temperature (day or night) used to make predictions.
#'
#' @export


preds_vs_obs_scatter_plot <- function(df, x, y, covar){

  x_range <- pretty(df[, x])
  y_range <- pretty(df[, y])

  lm <- lm(as.formula(paste0(y, "~ R0_1 - 1")), data = df)
  r_sq <- round(summary(lm)$r.squared, 3)

  par(mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0), xaxs = "r", yaxs = "r")

  plot(df[, x],
       df[, y],
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

  p <- recordPlot()

}
