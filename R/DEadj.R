# DEadj function: A helper function that shifts, scales, and rotates a 1D array.

DEadj <- function(par, par_hat, R2R1){
  par_adj <- (par - par_hat)%*%R2R1 + par_hat
  return(par_adj)
}
