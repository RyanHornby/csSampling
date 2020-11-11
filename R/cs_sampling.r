
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
  
  #Check weights
  #Check that the weights exist in both the survey object and the stan data
  if (is.null(weights(svydes))) {
    if (!is.null(weights(data_stan))) {
      warning("No survey weights, using weights from stan data instead")
      
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(weights(svydes))) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  #Check that the weights are the same
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), as.numeric(weights(svydes))))) {
    stop("Survey weights and stan data weights do not match")
  }
  #Check that the mean is 1
  if (mean(weights(data_stan)) != 1) {
    warning("Mean of the weights is not 1, scaling weights to a mean of one")
    new_weights = weights(data_stan)/mean(weights(data_stan))
    
    
  }
  
  
  print("stan fitting")
  out_stan  <- sampling(object = mod_stan, data = data_stan,
                        pars = par_stan,
                        chains = ctrl_stan$chains,
                        iter = ctrl_stan$iter, warmup = ctrl_stan$warmup, thin = ctrl_stan$thin
  )
  
  #Extract parameter draws and convert to unconstrained parameters
  
  #Get posterior mean (across all chains)
  par_samps_list <- rstan::extract(out_stan, pars = par_stan, permuted = TRUE)
  
  #concatenate across multiple chains
  par_samps <- as.matrix(out_stan, pars = par_stan)
  
  #convert to list type input > convert to unconstrained parameterization > back to matrix/array
  for(i in 1:dim(par_samps)[1]){
    if(i == 1){upar_samps <- unconstrain_pars(out_stan, list_2D_row_subset(par_samps_list, i))
    }else{upar_samps <- rbind(upar_samps, unconstrain_pars(out_stan, list_2D_row_subset(par_samps_list, i)))}
  }
  row.names(upar_samps) <- 1:dim(par_samps)[1]
  
  upar_hat <- colMeans(upar_samps) 
  
  #Estimate Hessian
  Hhat  <- -1*optimHess(upar_hat, gr = function(x){grad_log_prob(out_stan, x)})
  
  #create svrepdesign
  if(rep_design == TRUE){svyrep <- svydes
  }else{
    svyrep <- as.svrepdesign(design = svydes, type = ctrl_rep$type, replicates = ctrl_rep$replicates)
  }
  
  #Estimate Jhat = Var(gradient)
  print("gradient evaluation")
  rep_tmp <- withReplicates(design = svyrep, theta = grad_par, stanmod = mod_stan,
                            standata = data_stan, par_stan = par_stan, par_hat = upar_hat)#note upar_hat
  Jhat <- vcov(rep_tmp)
  
  #compute adjustment
  #use pivot for numerical stability - close to positive semi-definite if some parameters are highly correlated
  #(Q <- chol(m, pivot = TRUE))
  ## we can use this by
  #pivot <- attr(Q, "pivot")
  #Q[, order(pivot)]
  Hi <- solve(Hhat)
  V1 <- Hi%*%Jhat%*%Hi
  R1 <- chol(V1,pivot = TRUE)
  pivot <- attr(R1, "pivot")
  R1 <- R1[, order(pivot)]
  
  R2 <- chol(Hi, pivot = TRUE)
  pivot2 <- attr(R2, "pivot2")
  R2 <- R2[, order(pivot)]
  R2i <- solve(R2)
  R2iR1 <- R2i%*%R1
  
  #adjust samples
  upar_adj <- aaply(upar_samps, 1, DEadj, par_hat = upar_hat, R2R1 = R2iR1, .drop = TRUE)
  
  #back transform to constrained parameter space
  for(i in 1:dim(upar_adj)[1]){
    if(i == 1){par_adj <- unlist(constrain_pars(out_stan, upar_adj[i,])[par_stan])#drop derived quantities
    }else{par_adj <- rbind(par_adj, unlist(constrain_pars(out_stan, upar_adj[i,])[par_stan]))}
  }
  
  #make sure names are the same for sampled and adjusted parms
  row.names(par_adj) <- 1:dim(par_samps)[1]
  colnames(par_samps) <- colnames(par_adj)
  
  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))
  
  return(rtn)
  
}#end of cs_sampling

#' cs_sampling_brms
#'
#'
#' @import brms
#'
#'
#' @export
cs_sampling_brms <- function(svydes, mod, data, family, 
                             ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                             rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap")) {
  
  
  stancode <- make_stancode(mod, data = data, family = family)
  mod_brms  <- stan_model(model_code = stancode)
  data_brms <- make_standata(mod, data = data, family = family)
  par_brms<- c("b") #subset of parameters interested in
  
  return(cs_sampling(svydes = svy14rep, mod_stan = mod_brms, par_stan = par_brms, data_stan = data_brms, rep_design = TRUE))
  
}

#'
#' @import GGally
#' 
#' @method plot cs_sampling
#' @export
plot.cs_sampling <- function(x) {
  
  datpl <- data.frame(rbind(as.matrix(x$sampled_parms), as.matrix(x$adjusted_parms))
                      , as.factor(c(rep("NO", dim(x$sampled_parms)[1]), rep("YES", dim(x$adjusted_parms)[1]))))
  names(datpl)[dim(x$sampled_parms)[2]+1] <- c("Adjust")
  rownames(datpl) <- NULL
  
  
  require(GGally)
  
  my_ellipse <- function(data, mapping){
    ggplot(data = data, mapping = mapping) +
      geom_point()+
      stat_ellipse(level = 0.90, type = "norm", size = 2)
  }
  
  my_violin <- function(data, mapping){
    ggplot(data = data, mapping = mapping) +
      geom_violin(trim=TRUE,draw_quantiles = c(0.05, 0.5, 0.95),alpha=0.5, size = 1.5)
  }
  
  p1 <- ggpairs(datpl, mapping = aes(color = Adjust, alpha = 0.5), columns = c(1:(dim(datpl)[2]-1)),
                lower = list(continuous = my_ellipse))
  return(p1)
}


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
#' @import pkgcond
#' @export
grad_par <- function(pwts, svydata, stanmod, standata,par_stan,par_hat){
  #ignore svydata argument it allows access to svy object data
  standata$weights <- pwts
  
  
  suppress_messages(out_stan  <- sampling(object = stanmod, data = standata,
                        pars = par_stan,
                        chains = 0, warmup = 0,), "the number of chains is less than 1")
  
  gradpar <- grad_log_prob(out_stan,par_hat)
  return(gradpar)
}#end of grad theta


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

#helper function to convert list of MCMC output from stan to input for unconstrained tranformation


#' list_2D_row_subset
#'
#' @param nmlist - a named list whose elements are 1D or 2D arrays with the same number of rows
#' @param rindex - the row index to subset
#' @return a names list of types as the input, with elements subset to rows index provided
#' @export
list_2D_row_subset <- function(nmlist, rindex){
  temp_list <- list()
  k <- length(nmlist)
  for(k in 1:length(nmlist))
    if(is.na(dim(nmlist[[k]])[2])){#if the list element is only one dimension (not two)
      eval(parse(text = paste("temp_list$",names(nmlist)[k], " <- ", 
                              "(nmlist$",names(nmlist)[k], ")[rindex]", sep = "")))
    }else{eval(parse(text = paste("temp_list$",names(nmlist)[k], " <- ", 
                                  "(nmlist$",names(nmlist)[k], ")[rindex,]", sep = "")))
    }
  return(temp_list)
}
