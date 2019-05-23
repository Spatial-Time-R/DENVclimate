

#------------------------------------------------------------------------------
#' loop
#'
#' \code{loop} wrapper for \code{lapply}.
#'
#' @param parallel If TRUE it uses \code{parLapply}.

#' @export


loop <- function(..., parallel) {
  if (parallel) {
    parallel::parLapply(NULL, ...)
  } else {
    lapply(...)
  }
}


#------------------------------------------------------------------------------
#' loop_simplify
#'
#' \code{loop_simplify} simplify the output of \code{loop}.
#'
#' @param what the return value from \code{loop}.

#' @export


loop_simplify <- function(..., what) {
  vapply(loop(...), identity, what)
}


#------------------------------------------------------------------------------
#' write_out_rds
#'
#' \code{write_out_rds} saves an rds file.
#'
#' @param dat the dataframe to save.
#' @param my_path the output directory.
#' @param file_name the output file name.

#' @export


write_out_rds <- function(dat, my_path, file_name) {

  dir.create(my_path, FALSE, TRUE)

  saveRDS(dat, file.path(my_path, file_name))

}


#------------------------------------------------------------------------------
#' write_out_csv
#'
#' \code{write_out_csv} saves a csv file.
#'
#' @param dat the dataframe to save.
#' @param my_path the output directory.
#' @param file_name the output file name.

#' @export


write_out_csv <- function(dat, my_path, file_name) {

  dir.create(my_path, FALSE, TRUE)

  write.table(
    dat,
    file.path(my_path, file_name),
    row.names = FALSE,
    sep = ",")

}


#------------------------------------------------------------------------------
#' df_to_list
#'
#' \code{df_to_list} converts a dataframe into a list.
#'
#' @param x the dataframe to convert.
#' @param use_names if TRUE it saves the df column names into each list component.

#' @export


df_to_list <- function (x, use_names) {
  keep <- c("names", "class", "row.names")
  at <- attributes(x)
  attributes(x) <- at[intersect(names(at), keep)]
  ret <- unname(lapply(split(x, seq_len(nrow(x))), as.list))
  if (!use_names) {
    ret <- lapply(ret, unname)
  }
  if (is.character(at$row.names)) {
    names(ret) <- at$row.names
  }
  ret
}


#------------------------------------------------------------------------------
#' save_plot
#'
#' \code{save_plot} save a png file of a plot
#'
#' @param plot_obj The plot object
#' @param out_pth The path where to save the plot
#' @param out_fl_nm The output file name
#' @param wdt The plot width
#' @param hgt The plot height
#'
#' @export


save_plot <- function(plot_obj, out_pth, out_fl_nm, wdt, hgt){

  dir.create(out_pth, FALSE, TRUE)
  png(file.path(out_pth, paste0(out_fl_nm, ".png")),
      width = wdt,
      height = hgt,
      units = "cm",
      pointsize = 12,
      res = 300)
  print(plot_obj)
  on.exit(dev.off())

}


#------------------------------------------------------------------------------
#' basic_scatter_plot
#'
#' \code{basic_scatter_plot} save a png file of a plot
#'
#' @param df dataframe with the data to plot.
#' @param x variable to plot on the x axis.
#' @param y variable to plot on the y axis.
#'
#' @export


basic_scatter_plot <- function(df, x, y){

  par(mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))

  plot(df[, x],
       df[, y],
       xlab = x,
       ylab = y,
       pch = 19,
       cex = 0.5)

  p <- recordPlot()

}
