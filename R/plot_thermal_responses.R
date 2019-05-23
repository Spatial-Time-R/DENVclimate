

#------------------------------------------------------------------------------
#' plot_thermal_responses_all_data
#'
#' \code{plot_thermal_responses_all_data} makes a faceted plot of the
#'   thermal responses with all the observations.
#'
#' @param my_data Dataframe of observations.
#' @param all_fits Dataframe of model fits.
#'
#' @inheritParams save_plot
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon facet_wrap
#'   coord_cartesian scale_x_continuous
#'
#' @export


plot_thermal_responses_all_data <- function(my_data, all_fits, out_pth, out_fl_nm){

  p <- ggplot() +
    geom_point(aes(x = T, y = trait, colour = ref), data = my_data) +
    geom_line(aes(x = temp, y = mean), data = all_fits) +
    geom_ribbon(aes(x = temp, ymin = q1, ymax = q2), data = all_fits, alpha = 0.2) +
    facet_wrap(~ trait.name, scales = "free_y") +
    coord_cartesian(xlim = c(10, 40)) +
    scale_x_continuous("temperature (C)", breaks = seq(10,40,5), labels = seq(10,40,5))

  save_plot(p,
            out_pth,
            out_fl_nm = paste0(out_fl_nm, ".png"),
            wdt = 23,
            hgt = 15)

}


#------------------------------------------------------------------------------
#' plot_thermal_responses_mean_data
#'
#' \code{plot_thermal_responses_mean_data} makes a faceted plot of the
#'   thermal responses with the mean of the observations.
#'
#' @param my_data_av Dataframe of average observations.
#' @param all_fits Dataframe of model fits.
#'
#' @inheritParams save_plot
#'
#' @importFrom ggplot2 ggplot geom_pointrange aes geom_point geom_line geom_ribbon
#'  facet_wrap coord_cartesian scale_x_continuous scale_y_continuous
#'
#' @export


plot_thermal_responses_mean_data <- function(my_data_av, all_fits, out_pth, out_fl_nm){

  p <- ggplot() +
    geom_pointrange(aes(x = T, y = m, ymin = m - se, ymax = m + se), fatten = 2, size = 0.5, data = my_data_av) +
    geom_line(aes(x = temp, y = mean), data = all_fits) +
    geom_ribbon(aes(x = temp, ymin = q1, ymax = q2), data = all_fits, alpha = 0.2) +
    facet_wrap(~ trait.name, scales = "free_y") +
    coord_cartesian(xlim = c(10, 40)) +
    scale_y_continuous("trait") +
    scale_x_continuous("temperature (C)", breaks = seq(10,40,5), labels = seq(10,40,5))

  save_plot(p,
            out_pth,
            out_fl_nm = paste0(out_fl_nm, ".png"),
            wdt = 18,
            hgt = 15)

}
