
#------------------------------------------------------------------------------
#' myR0
#'
#' \code{myR0} estimates virus reproduction number R0 based on a Ross-MacDonald model.
#'
#' @param a mosquitoes biting rate.
#' @param b proportion of infectious bites that infect susceptible humans.
#' @param c proportion of bites on infected humans that infect previously
#'  uninfected mosquitoes.
#' @param PDR parasite development rate.
#' @param MDR mosquito immature development rate.
#' @param EFD number of eggs produced per female mosquito per day.
#' @param e2a mosquito egg-to-adult survival probability.
#' @param lf mosquito lifespan.
#'
#' @export


myR0 <- function(a, b, c, PDR, MDR, EFD, e2a, lf){

  # Creating a small constant to keep denominators from being zero.
  ec <- 0.000001

  mu <- 1 / (lf + ec)

  bc <- (b *c)

  ((a ^ 2 * bc * (EFD * e2a * MDR / (mu) ^ 2) * exp((-mu / (PDR + ec)))) / (mu)) ^ 0.5

}
