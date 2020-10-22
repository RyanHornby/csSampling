
#take in stan mod, list for stan data, name of parameters sampled,
#survey design object (or rep)
#return sampling object with parameters overwritten

#' cs_sampling
#'
#' cs_sampling is a wrapper function. It takes in a svydesign object (survey) and a stan model and inputs (rstan).
#' It calls the rstan::sampling to generate MCMC draws from the model. The constrained parameters are converted to unconstrained, adjusted, converted back to constrained and then output.
#' The adjustment process estimates a sandwich matrix adjustment to the posterior variance from two information matrices H and J. 
#' J is estimated via resampling with the withReplicates. For each set of replicate weights, rstan:sampling is called with no chains to instantiate a stanfit object. 
#' The stanfit object has an associated method grad_log_prob, which returns the gradient for a given input of unconstrained parameters. 
#' The variance of this gradient is taken across the replicates and provides a estimate of the J matrix
#' H is a Hessian estimated at the posterior mean via optimHess (stats) using grad_log_prob and the posterior mean as inputs
#' The asymptotic covariance for the posterior mean is then calculated as Hi V Hi, where Hi is the inverse of H. 
#' The asymptotic covariance for the posterior sampling procedure (due to mis-specification) is Hi. We take the "sqaure" root of these matrices
#' via cholesky decomposition R1'R1 = Hi V Hi and R2'R2 = Hi. 
#' The final adjustment rescales/rotates the posterior sample by R2iR1 where R2i is the inverse of R2. 
#' The final adjust can be interpreted as an asymptotic correction for model mis-specification due to using survey sampling weights as plug in values in the likelihood. This is often know as a "design effect" which is the "ratio" between the variance from simple random sample (Hi) and a complex survey sample (HiVHi)
#' @references Williams, M. R., and Savitsky, T. D. (2020) Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling. International Statistical Review, https://doi.org/10.1111/insr.12376. 
#'
#' @import rstan
#' @import survey
#' @import plyr
#' @param svydes - a svydesign object or a svyrepdesign object (see survey package). This contains cluster ID, strata, and weight information (svydesign) or replicate weight information (svyrepdesign)
#' @param mod_stan - a compiled stan model to be called by rstan::sampling
#' @param par_stan - a list of the parameter names in the parameter block of the stan model. These are treated as potentially constrained and will be converted to unconstrained parameters, adjusted, and then converted back  
#' @param data_stan - a list of data inputs for rstan::sampling() associated with mod_stan
#' @param ctrl_stan - a list of control parameters to pass to rstan::sampling(). Currently includes the number of chains, iter, warmpup, and thin with defualts
#' @param rep_design - logical indicating if the svydes object is a svyrepdesign. If FALSE, the design will be converted to a svyrepdesign using ctrl_rep settings
#' @param ctrl_rep - a list of settings when converting svydes from a svydesign object to a svyrepdesign object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate 
#' @return A list of the following:
#' \itemize{
#'  \item stan_fit - the original stanfit object returned by rstan::sampling for the weighted model
#'  \item sampled_parms - the array of parameters extracted from stan_fit corresponding to the parameter block in the stan model (specificed by stan_pars)
#'  \item adjusted_parms - the array (or matrix?) of adjusted parameters, corresponding to sampled_parms which have been rescaled and rotated
#' }
#' @export

cs_sampling <- function(svydes, mod_stan, par_stan, data_stan,
		ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
		rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap")){
  require(rstan)
  require(survey)
  require(plyr)
  
  #run STAN model
  print("stan fitting")
  out_stan  <- sampling(object = mod_stan, data = data_stan,
                            pars = par_stan,
                            chains = ctrl_stan$chains,
                            iter = ctrl_stan$iter, warmup = ctrl_stan$warmup, thin = ctrl_stan$thin
                            )
  
  #Get posterior mean (across all chains)
  par_samps <- extract(out_stan, pars = par_stan, permuted = FALSE)
  par_hat <- colMeans(par_samps, dim = 2)#dim = 1 by chain, dim = 2 across chains
  
  #Estimate Hessian
  Hhat  <- -1*optimHess(par_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  #create svrepdesign
  if(rep_design == TRUE){svyrep <- svydes
  	}else{
  	svyrep <- as.svrepdesign(design = svydes, type = ctrl_rep$type, replicates = ctrl_rep$replicates)
  }
  
  #Estimate Jhat = Var(gradient)
  print("gradient evaluation")
  rep_tmp <- withCallingHandlers(withReplicates(design = svyrep, theta = grad_par, stanmod = mod_stan,
  						standata = data_stan, par_stan = par_stan, par_hat = par_hat), warning = hideChainsWarnings)
  Jhat <- vcov(rep_tmp)
  
  #compute adjustment
  Hi <- solve(Hhat)
  V1 <- Hi%*%Jhat%*%Hi
  R1 <- chol(V1)
  R2i <- chol(Hi)
  R2 <- solve(R2i)
  R2R1 <- R2%*%R1
  
  #adjust samples
  par_adj <- aaply(par_samps, 1, DEadj, par_hat = par_hat, R2R1 = R2R1, .drop = FALSE)
  #matches par_samps if needed
  
  return(list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms =par_adj))

}#end of cs_sampling


##helper functions###

##grad_par helper function to nest within withReplicates()
#Stan will pass warnings from calling 0 chains, but still create out_stan object with
#grad_log_prob() method



#' grad_par
#'
#' @param par_hat - single set of unconstrained parameters (or estimates like posterior mean) to evaluate the gradient at  
#' @param pwts - is the weights argument. withReplicates updates this for each set of replicate weights
#' @param svydata - allows access to svrepdesign object's associated data. withReplicates expects it, but we do not use it. We use standata instead
#' @param stanmod - the compiled stan model to be passed to rstan::sampling
#' @param standata - list of data inputs to be passed to rstan::samlping
#' @param par_stan - list of potentially constrained parameters from the parameters block of the stan model
#' @return the gradient of the log posterior evaluated at par_hat
#' @export
grad_par <- function(pwts, svydata, stanmod, standata,par_stan,par_hat){
#ignore svydata argument it allows access to svy object data
standata$weights <- pwts




out_stan  <- sampling(object = stanmod, data = standata,
                      pars = par_stan,
                      chains = 0, warmup = 0,
                      )

gradpar <- grad_log_prob(out_stan,par_hat)
return(gradpar)
}#end of grad theta

hideChainsWarnings <- function(w) {
  if(any(grepl("the number of chains is less than 1", w))) 
    invokeRestart("muffleWarning")
}

#helper function to apply matrix rotation

#' DEadj
#'
#' @param par - a 1D array/matrix of parameter values to be rescaled/rotated
#' @param par_hat - a 1D array/matrix of center values, typically the mean
#' @param R2R1 - the rotation matrix to rescale/rotate par
#' @return A 1D matrix corresponding to the input par that has been rescaled/rotated
#' @export
DEadj <- function(par, par_hat, R2R1){
par_adj <- (par - par_hat)%*%R2R1 + par_hat
return(par_adj)
}