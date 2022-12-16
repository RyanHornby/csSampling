# grad_par: A helper function to use with survey::withReplicates, which estimates the gradient of the log of Stan probability model for a given set of weights.
#' @import pkgcond
#'
grad_par <- function(pwts, svydata, stanmod, standata,par_hat){
  #ignore svydata argument - it allows access to svy object data
  standata$weights <- pwts


  pkgcond::suppress_messages(out_stan  <- rstan::sampling(object = stanmod, data = standata,
                                          chains = 0, warmup = 0,), "the number of chains is less than 1")

  gradpar <- rstan::grad_log_prob(out_stan,par_hat)
  return(gradpar)
}
