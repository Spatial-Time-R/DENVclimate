#------------------------------------------------
#' mean_trait
#'
#' \code{mean_trait} calculates the mean
#'  of the trait values and standard deviation
#'
#' @param my_data Dataframe of observations.
#'
#' @importFrom dplyr group_by %>% summarise n
#'
#' @export


mean_trait <- function(my_data){

  my_data %>%
  group_by(trait.name, T) %>%
  summarise(m = mean(trait), sd = sd(trait), n = n())

}
