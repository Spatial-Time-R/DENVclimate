get_covariate_scaling_factor <- function(var, ppyear, year.f, year.i){
  
  scale <- 1
  
  if(grepl("Re.", var) | grepl("Im.", var)){
    
    scale <- ppyear * (year.f - year.i + 1) / 2 
    
  } 
  
  if(grepl("const_term$", var)){
    
    scale <- ppyear * (year.f - year.i + 1) 
    
  }  
  
  scale
  
}
