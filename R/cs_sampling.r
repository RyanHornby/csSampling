
#take in stan mod, list for stan data, name of parameters sampled,
#survey design object (or rep)
#return sampling object with parameters overwritten

#' cs_sampling
#'
#' cs_sampling is a wrapper function. It takes in a \code{\link[survey]{svydesign}} object and a \code{\link[rstan]{stan_model}} and inputs for \code{\link[rstan]{sampling}}.
#' It calls the \code{\link[rstan]{sampling}} to generate MCMC draws from the model. The constrained parameters are converted to unconstrained, adjusted, converted back to constrained and then output.
#' The adjustment process estimates a sandwich matrix adjustment to the posterior variance from two information matrices H and J. 
#' J is estimated via resampling with \code{\link[survey]{withReplicates}}. For each set of replicate weights, \code{\link[rstan]{sampling}} is called with no chains to instantiate a \code{\link[rstan]{stanfit-class}} object. 
#' The \code{\link[rstan]{stanfit-class}} object has an associated method \code{\link[rstan]{grad_log_prob}}, which returns the gradient for a given input of unconstrained parameters. 
#' The variance of this gradient is taken across the replicates and provides a estimate of the J matrix.
#' H is a Hessian estimated at the posterior mean via \code{\link[stats]{optimHess}} using \code{\link[rstan]{grad_log_prob}} and the posterior mean as inputs.
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
#' @import pkgcond
#' @param svydes - a \code{\link[survey]{svydesign}} object or a \code{\link[survey]{svrepdesign}} object. This contains cluster ID, strata, and weight information (\code{\link[survey]{svydesign}}) or replicate weight information (\code{\link[survey]{svrepdesign}})
#' @param mod_stan - a compiled stan model to be called by \code{\link[rstan]{sampling}}
#' @param par_stan - a list of a subset of parameters to output after adjustment. All parameters are adjusted including the derived parameters, so users may want to only compare subsets. The default, NA, will return all parameters.  
#' @param data_stan - a list of data inputs for \code{\link[rstan]{sampling}} associated with mod_stan
#' @param ctrl_stan - a list of control parameters to pass to \code{\link[rstan]{sampling}}. Currently includes the number of chains, iter, warmup, and thin with defaults
#' @param rep_design - logical indicating if the svydes object is a \code{\link[survey]{svrepdesign}}. If FALSE, the design will be converted to a \code{\link[survey]{svrepdesign}} using ctrl_rep settings
#' @param ctrl_rep - a list of settings when converting svydes from a \code{\link[survey]{svydesign}} object to a \code{\link[survey]{svrepdesign}} object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate (see \code{\link[survey]{as.svrepdesign}}). 
#' @param sampling_args - a list of extra arguments that get passed to \code{\link[rstan]{sampling}}.
#' @return A list of the following:
#' \itemize{
#'  \item stan_fit - the original \code{\link[rstan]{stanfit-class}} object returned by \code{\link[rstan]{sampling}} for the weighted model
#'  \item sampled_parms - the array of parameters extracted from stan_fit corresponding to the parameter block in the stan model (specified by stan_pars)
#'  \item adjusted_parms - the array of adjusted parameters, corresponding to sampled_parms which have been rescaled and rotated.
#' }
#' @examples 
#'###Custom Stan Model###
#'#Weighted dirichlet-multinomial model
#'
#'#survey package example
#'data(api)
#'##make sure weights sum to n##
#'apiclus1$newpw <- apiclus1$pw/mean(apiclus1$pw)
#'
#'dclus1<-svydesign(id=~dnum, weights=~newpw, data=apiclus1, fpc=~fpc)
#'svymean(~stype, dclus1)
#'
#'#Use default replicate design
#'rclus1<-as.svrepdesign(dclus1)
#'svymean(~stype, rclus1)
#'
#'#use cs_sampling
#'mod_dm <- stan_model(model_code = csSampling::proportion_estimate)
#'
#'#Set the Data for Stan
#'y <- as.factor(rclus1$variables$stype)
#'yM <- model.matrix(~y -1)
#'n <- dim(yM)[1]
#'K <- dim(yM)[2]
#'#Uniform prior for alpha for now
#'alpha<-rep(1,K)
#'weights <- rclus1$pweights
#'
#'#Create stan Data list
#'data_stan<-list("y"=yM,"alpha"=alpha,"K"=K, "n" = n, "weights" = weights)
#'ctrl_stan<-list("chains"=1,"iter"=2000,"warmup"=1000,"thin"=1)
#'
#'mod1 <- cs_sampling(svydes = rclus1, mod_stan = mod_dm, data_stan = data_stan, ctrl_stan = ctrl_stan, rep_design = T)
#'
#'#plot parameters of interest - proportions (thetas)
#'plot(mod1, varnames = paste("theta",1:3, sep = ""))
#'#compare summary statistics
#'#unadjusted weighted Bayes
#'cbind(colMeans(mod1$sampled_parms[,4:6]), sqrt(diag(cov(mod1$sampled_parms[,4:6]))))
#'#hybrid variance adjusted
#'cbind(colMeans(mod1$adjusted_parms[,4:6]), sqrt(diag(cov(mod1$adjusted_parms[,4:6]))))
#'#Taylor linearization
#'svymean(~stype, dclus1)
#'#Replication
#'svymean(~stype, rclus1)
#' @export

cs_sampling <- function(svydes, mod_stan, par_stan = NA, data_stan,
                        ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                        rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"), 
                        sampling_args = list()){
  require(rstan)
  require(survey)
  require(plyr)
  require(pkgcond)
  
  #Check weights
  #Check that the weights exist in both the survey object and the stan data
  #weights() returns full replicate weights set if svrepdesign
  if(rep_design){svyweights <- svydes$pweights}else{svyweights <-weights(svydes)}
  
  if (is.null(svyweights)) {
    if (!is.null(weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  #Check that the weights are the same
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  #Check that the mean is 1
  if (mean(weights(data_stan)) != 1) {
    stop("Mean of the weights is not 1")
  }
  
  
  print("stan fitting")
  out_stan  <- do.call(sampling, c(list(object = mod_stan, data = data_stan,
                                        pars = par_stan,
                                        chains = ctrl_stan$chains,
                                        iter = ctrl_stan$iter, warmup = ctrl_stan$warmup, thin = ctrl_stan$thin), sampling_args)
  )
  
  #Extract parameter draws and convert to unconstrained parameters
  
  #Get posterior mean (across all chains)
  par_samps_list <- rstan::extract(out_stan, permuted = TRUE)
  
  #If par_stan is not provided (NA) use all parameters (except "lp__", which is last)
  if(anyNA(par_stan)){
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }
  
  #concatenate across multiple chains - save for later for export
  par_samps <- as.matrix(out_stan, pars = par_stan)
  
  #convert to list type input > convert to unconstrained parameterization > back to matrix/array
  for(i in 1:dim(par_samps)[1]){#just need the length here
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
                            standata = data_stan, par_hat = upar_hat)#note upar_hat
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
  pivot2 <- attr(R2, "pivot")
  R2 <- R2[, order(pivot2)]
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
#'cs_sampling_brms is a wrapper function that takes inputs in the form of model statements in familiar brms syntax.
#'Then brms helper functions build Stan models and call cs_sampling.
#'
#'
#' @import brms
#' @param svydes - a \code{\link[survey]{svydesign}} object or a \code{\link[survey]{svrepdesign}} object. This contains cluster ID, strata, and weight information (\code{\link[survey]{svydesign}}) or replicate weight information (\code{\link[survey]{svrepdesign}})
#' @param brmsmod - \code{\link[brms]{brmsformula}}  object, as input to \code{\link[brms]{make_stancode}}. The \code{\link[brms]{brmsformula}}  must specify a weight variable via weights().
#' @param par_brms - a list of a subset of parameters to output after adjustment. All parameters are adjusted including the derived parameters, so users may want to only compare subsets. The default, NA, will return all parameters.  
#' @param data - a data frame, as input to \code{\link[brms]{make_stancode}}
#' @param family - \code{\link[stats]{brmsfamily}} or \code{\link[brms]{brmsfamily}} as input to \code{\link[brms]{make_stancode}} specifying distribution and link function
#' @param prior - optional input to \code{\link[brms]{make_stancode}}
#' @param stanvars - optional input to \code{\link[brms]{make_stancode}}
#' @param knots - optional input to \code{\link[brms]{make_stancode}}
#' @param ctrl_stan - a list of control parameters to pass to \code{\link[rstan]{sampling}}. Currently includes the number of chains, iter, warmpup, and thin with defaults.
#' @param rep_design - logical indicating if the svydes object is a \code{\link[survey]{svrepdesign}}. If FALSE, the design will be converted to a \code{\link[survey]{svrepdesign}} using ctrl_rep settings
#' @param ctrl_rep - a list of settings when converting svydes from a \code{\link[survey]{svydesign}} object to a \code{\link[survey]{svrepdesign}} object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate (see \code{\link[survey]{as.svrepdesign}}). 
#' @param stancode_args - a list of extra arguments to be passed to \code{\link[brms]{make_stancode}}.
#' @param standata_args - a list of extra arguments to be passed to \code{\link[brms]{make_standata}}.
#' @param sampling_args - a list of extra arguments to be passed to \code{\link[rstan]{sampling}}.
#' @return The output of cs_sampling.
#' @examples
#' 
#' ####BRMS Wrapper#####
#' #Linear regression from survey package api data
#' data(api)
#' dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
#' ##need to scale the weights in the survey design and in the stan model##
#' apistrat$wtsscl <- apistrat$pw *length(apistrat$pw)/sum(apistrat$pw)
#'
#' dstrat_sc<-svydesign(id=~1,strata=~stype, weights=~wtsscl, data=apistrat, fpc=~fpc)
#' ###example 1 api00~ell+meals+mobility##
#' #Use cs_sampling directly#  
#' 
#'stancode <- make_stancode(brmsformula(api00|weights(wtsscl) ~ ell+meals+mobility,center = FALSE), 
#'                          data = apistrat, family = gaussian(), save_model = "brms_wt_lm.stan")
#'
#'mod_brms  <- stan_model('brms_wt_lm.stan')
#'
#'data_brms <- make_standata(brmsformula(api00|weights(wtsscl) ~ ell+meals+mobility,center = FALSE), 
#'                           data = apistrat, family = gaussian())
#'             
#'set.seed(12345)
#'mod.brms_1 <- cs_sampling(svydes = dstrat_sc, mod_stan = mod_brms, data_stan = data_brms, 
#'                          ctrl_stan = list(chains = 2, iter = 2000, warmup = 1000, thin = 2))
#'
#' #Compare to Wrapper#                          
#'set.seed(12345)
#'mod.brms_2 <- cs_sampling_brms(svydes = dstrat_sc, 
#'                        brmsmod = brmsformula(api00|weights(wtsscl) ~ ell+meals+mobility,center = FALSE), 
#'                        data = apistrat, family = gaussian())
#'                     
#' #compare to svyglm
#'summary(svyglm(api00~ell+meals+mobility, design=dstrat_sc))         
#'
#' #plot all parameters by default
#'plot(mod.brms_3)
#subset plot by varnames
#'plot(mod.brms_3, varnames = paste("b", 1:3, sep =""))
#'
#'pp <- plot(mod.brms_3)
#'pp[2,1]
#' @export
cs_sampling_brms <- function(svydes, brmsmod, data, family, par_brms = NA,prior = NULL, stanvars = NULL, knots = NULL, 
                             ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                             rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                             stancode_args = list(), standata_args = list(), sampling_args = list()) {
  
  
  stancode <- do.call(make_stancode, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), stancode_args))
  print("compiling stan model")
  mod_brms  <- stan_model(model_code = stancode)
  data_brms <- do.call(make_standata, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), standata_args))
  
  return(cs_sampling(svydes = svydes, mod_stan = mod_brms, par_stan = par_brms, data_stan = data_brms, 
                     rep_design = rep_design, ctrl_rep = ctrl_rep, ctrl_stan = ctrl_stan, sampling_args))
  
}


#' plot.cs_sampling
#'
#' Produce pairwise plots comparing each parameter distribution before and after adjustment. 
#' Takes the output list from cs_sampling as input.
#'
#' @import GGally
#' @param x - object of type cs_sampling
#' @param varnames - optional vector of names of subset of variable for pairs plotting
#' @return The output of \code{\link[GGally]{ggpairs}}.
#' 
#' @method plot cs_sampling
#' @export
plot.cs_sampling <- function(x, varnames = NULL) {
  
  datpl <- data.frame(rbind(as.matrix(x$sampled_parms), as.matrix(x$adjusted_parms))
                      , as.factor(c(rep("NO", dim(x$sampled_parms)[1]), rep("YES", dim(x$adjusted_parms)[1]))))
  names(datpl)[dim(x$sampled_parms)[2]+1] <- c("Adjust")
  rownames(datpl) <- NULL
  
  #subset to varnames
  if(!is.null(varnames)){datpl <- datpl[, c(varnames, "Adjust")]}
  
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
#' A helper function to use with \code{\link[survey]{withReplicates}}, which estimates the gradient of the log of Stan probability model for a given set of weights.
#'
#' @param par_hat - single set of unconstrained parameters (or estimates like posterior mean) to evaluate the gradient at  
#' @param pwts - the weights argument. \code{\link[survey]{withReplicates}} updates this for each set of replicate weights.
#' @param svydata - allows access to \code{\link[survey]{svrepdesign}} object's associated data. \code{\link[survey]{withReplicates}} expects it, but we do not use it. We use standata instead
#' @param stanmod - the compiled stan model to be passed to \code{\link[rstan]{sampling}}
#' @param standata - list of data inputs to be passed to \code{\link[rstan]{sampling}}
#' @return the gradient of the log posterior for stanmod evaluated at par_hat
#' @import pkgcond
#' @export
grad_par <- function(pwts, svydata, stanmod, standata,par_hat){
  #ignore svydata argument - it allows access to svy object data
  standata$weights <- pwts
  
  
  suppress_messages(out_stan  <- sampling(object = stanmod, data = standata,
                                          chains = 0, warmup = 0,), "the number of chains is less than 1")
  
  gradpar <- grad_log_prob(out_stan,par_hat)
  return(gradpar)
}#end of grad theta


#helper function to apply matrix rotation

#' DEadj
#' 
#' A helper function that shifts, scales, and rotates a 1D array.
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

#helper function to convert list of MCMC output from stan to input for unconstrained transformation


#' list_2D_row_subset
#' 
#' A helper function to convert a list of MCMC output from Stan to input comparable with the unconstrained transformation.
#'
#' @param nmlist - a named list whose elements are 1D or 2D arrays with the same number of rows
#' @param rindex - the row index to subset
#' @return a named list of types as the input, with elements subset to rows index provided
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
