myR0<-function(a, b, c, PDR, MDR, EFD, e2a, lf){
  ec <- 0.000001
  
  mu = 1/(lf + ec)
  bc = (b*c)
  ((a^2*bc*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5
}
