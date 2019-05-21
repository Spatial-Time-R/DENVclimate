loop <- function(..., parallel) {
  if (parallel) {
    parallel::parLapply(NULL, ...)
  } else {
    lapply(...)
  }
}

loop_simplify <- function(..., what) {
  vapply(loop(...), identity, what)
}

write_out_rds <- function(dat, my_path, file_name) {
  
  dir.create(my_path, FALSE, TRUE)
  
  saveRDS(dat, file.path(my_path, file_name))
  
}

write_out_csv <- function(dat, my_path, file_name) {
  
  dir.create(my_path, FALSE, TRUE)
  
  write.table(
    dat, 
    file.path(my_path, file_name), 
    row.names = FALSE, 
    sep = ",")
  
}

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

lm_eqn <- function(df, y, x){
  
  # fit linear model - no intercept 
  frmla <- as.formula(paste0(y, " ~ ", x, " - 1"))
  m <- lm(frmla, df)
  
  #browser()
  
  # calculates R2
  r2 <- 1 - crossprod(residuals(m)) / crossprod(df[,y] - mean(df[,y]))
  #r2 <- summary(m)$r.squared # R 
  
  eq <- substitute(italic(y) == b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(b = format(coef(m), digits = 4), 
                        r2 = format(r2, digits = 4)))
  as.character(as.expression(eq))                 
}
