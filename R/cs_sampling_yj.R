#' cs_sampling_yj
#'
#' \code{cs_sampling} is a wrapper function. It takes in a \code{\link[survey]{svydesign}} object and a \code{\link[rstan]{stan_model}} and inputs for \code{\link[rstan]{sampling}}.
#' It calls the \code{\link[rstan]{sampling}} to generate MCMC draws from the model. The constrained parameters are converted to unconstrained, adjusted, converted back to constrained and then output.
#' The adjustment process estimates a sandwich matrix adjustment to the posterior variance from two information matrices H and J.
#' J is estimated via resampling with \code{\link[survey]{withReplicates}}. For each set of replicate weights, \code{\link[rstan]{sampling}} is called with no chains to instantiate a \code{\link[rstan]{stanfit-class}} object.
#' The \code{\link[rstan]{stanfit-class}} object has an associated method \code{\link[rstan]{grad_log_prob}}, which returns the gradient for a given input of unconstrained parameters.
#' The variance of this gradient is taken across the replicates and provides a estimate of the J matrix.
#' '\code{cs_sampling} allows for different options for estimation of the Hessian matrix H.
#' By default H is estimated as the Monte Carlo mean via \code{\link[stats]{optimHess}} using \code{\link[rstan]{grad_log_prob}} and each posterior draw as inputs. Using just the posterior mean is faster but less stable (previous default).
#' The asymptotic covariance for the posterior mean is then calculated as Hi V Hi, where Hi is the inverse of H.
#' The asymptotic covariance for the posterior sampling procedure (due to mis-specification) is Hi.
#' By default, \code{cs_sampling} takes the "square" root of these matrices
#' via eigenvalue decomposition R1'R1 = Hi V Hi and R2'R2 = Hi. This is more stable but slower than using the Cholesky decomposition (previous default).
#' The final adjustment rescales/rotates the posterior sample by R2iR1 where R2i is the inverse of R2.
#' The final adjust can be interpreted as an asymptotic correction for model mis-specification due to using survey sampling weights as plug in values in the likelihood. This is often know as a "design effect" which is the "ratio" between the variance from simple random sample (Hi) and a complex survey sample (HiVHi).
#'
#' @references Williams, M. R., and Savitsky, T. D. (2020) Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling. International Statistical Review, https://doi.org/10.1111/insr.12376.
#'
#' @author Matt Williams.
#'
#' @param svydes - a \code{\link[survey]{svydesign}} object or a \code{\link[survey]{svrepdesign}} object. This contains cluster ID, strata, and weight information (\code{\link[survey]{svydesign}}) or replicate weight information (\code{\link[survey]{svrepdesign}})
#'
#' @param mod_stan - a compiled stan model to be called by \code{\link[rstan]{sampling}}
#'
#' @param par_stan - a list of a subset of parameters to output after adjustment. All parameters are adjusted including the derived parameters, so users may want to only compare subsets. The default, NA, will return all parameters.
#'
#' @param data_stan - a list of data inputs for \code{\link[rstan]{sampling}} associated with mod_stan
#'
#' @param ctrl_stan - a list of control parameters to pass to \code{\link[rstan]{sampling}}. Currently includes the number of chains, iter, warmup, and thin with defaults
#'
#' @param rep_design - logical indicating if the svydes object is a \code{\link[survey]{svrepdesign}}. If FALSE, the design will be converted to a \code{\link[survey]{svrepdesign}} using ctrl_rep settings
#'
#' @param ctrl_rep - a list of settings when converting svydes from a \code{\link[survey]{svydesign}} object to a \code{\link[survey]{svrepdesign}} object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate (see \code{\link[survey]{as.svrepdesign}}).
#'
#' @param matrix_sqrt - a string indicating the method to use to take the "square root" of the R1 and R2 matrices. The default "eigen" uses the eigenvalue decomposition. Otherwise, the Cholesky decomposition is used.
#'
#' @param diag_only - a logical indicating whether the variance adjustment should only use diagonals of H and J.
#'
#' @param subset_matrix - an index of (unconstrained) parameters to subset the adjustment. Requires knowledge of the stan model parameterization.
#' 
#' @param prior_only - a logical indicating if the stan model has an option for sampling just from the prior distribution. This can be used to further refine the estimates for covariances H and J.
#'
#' @param export_unconst_pars - a logical indicating whether to return the unconstrained parameters, both original and adjusted.
#' 
#' @param yj_range - a vector of the form c(lower, upper) to provide a range for the Yeo-Johnson transformation for parameter usually between (-3,3). Narrower values provide stability for inverting the transformation on the adjusted parameters.
#'
#' @param sampling_args - a list of extra arguments that get passed to \code{\link[rstan]{sampling}}.
#'
#' @import rstan
#' @import survey
#' @import plyr
#' @import pkgcond
#'
#'
#' 
#'
#'
#' @return A list of the following:
#' \itemize{
#'  \item stan_fit - the original \code{\link[rstan]{stanfit-class}} object returned by \code{\link[rstan]{sampling}} for the weighted model
#'  \item sampled_parms - the array of parameters extracted from stan_fit corresponding to the parameter block in the stan model (specified by stan_pars)
#'  \item adjusted_parms - the array of adjusted parameters, corresponding to sampled_parms which have been rescaled and rotated.
#' }
#'
#' @export
cs_sampling_yj <- function(svydes, mod_stan, par_stan = NA, data_stan,
                        ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                        rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                        matrix_sqrt = "eigen",
                        diag_only = FALSE,
                        subset_matrix = NULL,
                        prior_only = FALSE,
                        export_unconst_pars = FALSE,
                        yj_range = c(-2.5,2.5),
                        sampling_args = list()){

  .cs_sampling_yj_process(
    svydes = svydes,
    mod_stan = mod_stan,
    par_stan = par_stan,
    data_stan = data_stan,
    ctrl_stan = ctrl_stan,
    rep_design = rep_design,
    ctrl_rep = ctrl_rep,
    matrix_sqrt = matrix_sqrt,
    diag_only = diag_only,
    subset_matrix = subset_matrix,
    prior_only = prior_only,
    export_unconst_pars = export_unconst_pars,
    yj_range = yj_range,
    sampling_args = sampling_args
  )
}

.cs_sampling_yj_process <- function(svydes, mod_stan, par_stan = NA, data_stan,
                                    ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                                    rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                                    matrix_sqrt = "eigen",
                                    diag_only = FALSE,
                                    subset_matrix = NULL,
                                    prior_only = FALSE,
                                    export_unconst_pars = FALSE,
                                    yj_range = c(-2.5,2.5),
                                    sampling_args = list(),
                                    stan_fit = NULL){

  #Check weights
  #Check that the weights exist in both the survey object and the stan data
  #weights() returns full replicate weights set if svrepdesign
  if(rep_design){svyweights <- svydes$pweights}else{svyweights <-stats::weights(svydes)}

  if (is.null(svyweights)) {
    if (!is.null(stats::weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(stats::weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = stats::weights(svydes)
    }
  }
  #Check that the weights are the same
  if (!isTRUE(all.equal(as.numeric(stats::weights(data_stan)), as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  #Check that weights sum to the sample size
  if (abs(sum(stats::weights(data_stan)) -  length(stats::weights(data_stan))) > 1.0) {
    warning("Sum of the weights may not equal the sample size")
  }

  
  #Estimate Hessian for Prior
  if(prior_only){
  print("(0) Setting up Prior-Only Model (0)")
    prior_data <- data_stan
    prior_data$prior_only <- 1 #should convert to integer of 1
    pkgcond::suppress_messages(out_stan_prior  <- rstan::sampling(object = mod_stan, data = prior_data,
                                                                  chains = 0, warmup = 0,), "the number of chains is less than 1")
    
  }

  if (is.null(stan_fit)) {
    print("(1) stan fitting (1)")
    out_stan  <- do.call(rstan::sampling, c(list(object = mod_stan, data = data_stan,
                                                 pars = par_stan,
                                                 chains = ctrl_stan$chains,
                                                 iter = ctrl_stan$iter, warmup = ctrl_stan$warmup, thin = ctrl_stan$thin), sampling_args)
    )
  } else {
    print("(1) loading fitted stan model (1)")
    out_stan <- stan_fit
  }

  #Extract parameter draws and convert to unconstrained parameters

  #Get posterior mean (across all chains)
  par_samps_list <- rstan::extract(out_stan, permuted = TRUE)

  #If par_stan is not provided (NA) use all parameters (except "lp__", which is last)
  if(anyNA(par_stan)){
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }

  #concatenate across multiple chains - save for later for export
  par_samps <- as.matrix(out_stan, pars = par_stan)

  #number of MCMC draws
  ndraws <- dim(par_samps)[1]

  print("(2) Transforming Parameters (2)")
  #convert to list type input > convert to unconstrained parameterization > back to matrix/array
  ##This is the bottleneck###

  #preallocate upar_samps ahead of time instead of using rbind
  tmplist <- list_2D_row_subset(par_samps_list, 1)
  upar_samps_init <- unconstrain_pars(out_stan, tmplist)
  upar_samps <- matrix(data = NA, nrow = ndraws, ncol = length(upar_samps_init))

#transform parameters
 
    for(i in 1:ndraws){#just need the length here
      if(i %% 500 == 0){print(paste0("Converting draw ", i))} #status message for user
      tmplist <- list_2D_row_subset(par_samps_list, i)
      upar_samps[i,] <- rstan::unconstrain_pars(out_stan, tmplist)
    }
  
  
  #convert using yj transform - iterate across variables
  yju <- upar_samps
  lamvec <- rep(1, dim(upar_samps)[2])
  for(k in 1:(dim(upar_samps)[2])){
    tmpPower <- car::powerTransform(object = upar_samps[,k], family = "yjPower")
    lamvec[k] <- max(yj_range[1],min(yj_range[2],stats::coef(tmpPower, round = TRUE))) #make sure between -3 and 3
    yju[,k] <- VGAM::yeo.johnson(upar_samps[,k], lamvec[k])
  }
  
  
  #functions to convert
  der_yjinv <- function(x, lambda){numDeriv::grad(func = VGAM::yeo.johnson, x = x, lambda = lambda, inverse =TRUE)}
  
  grad_yjinv <- function(yj, lambda){#one MCMC draw at a time
    gradtmp <- rep(0, length(yj))
    for(k in 1:length(yj)){
      gradtmp[k] <- der_yjinv(x = yj[k], lambda = lambda[k])
    }
    return(gradtmp)
  }

  row.names(upar_samps) <- 1:dim(par_samps)[1]
  
  #posterior mean on transformed scale, then transform back
  yju_hat <- colMeans(yju)
  PV_yj <- stats::var(yju)
  
  Hhat <- solve(PV_yj)
  
  upar_hat <- yju_hat
  
  for(k in 1:(dim(upar_samps)[2])){
    upar_hat[k] <- VGAM::yeo.johnson(yju_hat[k], lamvec[k], inverse =TRUE)
  }


  #Estimate Hessian for Prior
  H0 <- NULL
  if(prior_only){ #we could also take the MCMC average but start simple here.
  H0  <- -1*stats::optimHess(upar_hat, gr = function(x){rstan::grad_log_prob(out_stan_prior, x)})
  }
  
  #create svrepdesign
  if(rep_design == TRUE){svyrep <- svydes
  }else{
    svyrep <- survey::as.svrepdesign(design = svydes, type = ctrl_rep$type, replicates = ctrl_rep$replicates)
  }

  #Estimate Jhat = Var(gradient)
  print("(3) Estimating Replicate Variance (3)")
  #perhaps a slowdown for large number of samples/data?
  rep_tmp <- survey::withReplicates(design = svyrep, theta = grad_par, stanmod = mod_stan,
                                    standata = data_stan, par_hat = upar_hat)#note upar_hat
  
  gvtmp <- grad_yjinv(yju_hat, lamvec)
  GMat <- t(t(gvtmp))%*%t(gvtmp)
  
  Jhat <- GMat*stats::vcov(rep_tmp)
  
  if(prior_only){ #non-asymptotic correction for prior
    Jhat <- Jhat + GMat*H0
  }

  print("(4) Estimating Adjustment (4)")
  #compute adjustment

  #independent or simultaneous adjustment
  if(diag_only){
  Hhat <- Matrix::Diagonal(n = dim(Hhat)[1], x = diag(Hhat))
  Jhat <- Matrix::Diagonal(n = dim(Jhat)[1], x = diag(Jhat)) 
  }
  #only adjust a subset of the unconstrained parameters - requires specific knowledge of stan model
  if(!is.null(subset_matrix)){#subset_matrix is an index of parameters (e.g. global)
    Htmp <- Hhat[subset_matrix, subset_matrix]
    Jtmp <- Jhat[subset_matrix, subset_matrix]
    ktmp <- dim(Hhat)[1]
    k1tmp <- dim(Htmp)[1]
    k2tmp <- ktmp - k1tmp
    Itmp <- Matrix::Diagonal(n = k2tmp, 1)
    Ztmp <- Matrix::Matrix(0, nrow = k1tmp, ncol = k2tmp, sparse = TRUE)
    tZtmp <- Matrix::Matrix(0, nrow = k2tmp, ncol = k1tmp, sparse = TRUE)
    
    Hhat <- rbind(
              cbind(Htmp, Ztmp),
              cbind(tZtmp, Itmp)
    )
    
    Jhat <- rbind(
      cbind(Jtmp, Ztmp),
      cbind(tZtmp, Itmp)
    )
    
  }
  
  Hi <- solve(Hhat)
  V1 <- Hi%*%Jhat%*%Hi

  if(matrix_sqrt == "eigen"){#use eigenvalue decomposition
    eigV <- eigen(V1, symmetric = TRUE)
    R1 <- diag(sqrt(abs(eigV$values)))%*%t(eigV$vectors)

    eigHi <- eigen(Hi, symmetric = TRUE)
    R2 <- diag(sqrt(abs(eigHi$values)))%*%t(eigHi$vectors)
  }else{#use cholesky decomposition
    R1 <- chol(V1,pivot = TRUE)
    pivot <- attr(R1, "pivot")
    R1 <- R1[, order(pivot)]

    R2 <- chol(Hi, pivot = TRUE)
    pivot2 <- attr(R2, "pivot")
    R2 <- R2[, order(pivot2)]
  }
  R2i <- solve(R2)
  R2iR1 <- R2i%*%R1

  #adjust samples
  print("(5) Applying Adjustment (5)")
  yju_adj <- plyr::aaply(yju, 1, DEadj, par_hat = yju_hat, R2R1 = R2iR1, .drop = TRUE)

  #back transform to constrained parameter space
  #special cases (4) needed for combinations of constrained and unconstrained par have 1 dimension
  #not likely/possible for unconstrained dim > constrained so one case might never be used
  

  print("(6) Backtransforming Parameters (6)")
  
  #invert yj transform
  upar_adj <- yju_adj
  
  for(k in 1:(dim(upar_samps)[2])){
    upar_adj[,k] <- VGAM::yeo.johnson(yju_adj[,k], lamvec[k], inverse = TRUE)
  }
  
  #convert back  to constrained parameters
  if(is.null(dim(upar_adj))){upardim <- length(upar_adj)
    upardim <- dim(upar_adj)[1] #treat 1 dimensional parameter as special due to dimension drop
    par_adj_tmp <- unlist(rstan::constrain_pars(out_stan, upar_adj[1])[par_stan])
    par_adj <- matrix(data = NA, nrow = upardim, ncol = length(par_adj_tmp))
    colnames(par_adj) <- names(par_adj_tmp)
    for (i in 1:upardim) {
      if(i %% 500 == 0){print(paste0("Back Converting draw ", i))} #status message for user
      if(upardim == 1){par_adj[i] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i])[par_stan])
      }else{par_adj[i,] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i])[par_stan])}
    }#treat 1 dimensional parameter as special due to dimension drop
  }else{
    upardim <- dim(upar_adj)[1] #treat 1 dimensional parameter as special due to dimension drop
    #only difference between top and bottom is the comma [i,] for different dimensions
    par_adj_tmp <- unlist(rstan::constrain_pars(out_stan, upar_adj[1,])[par_stan])
    par_adj <- matrix(data = NA, nrow = upardim, ncol = length(par_adj_tmp))
    colnames(par_adj) <- names(par_adj_tmp)
    for (i in 1:upardim) {
      if(i %% 500 == 0){print(paste0("Back Converting draw ", i))} #status message for user
      if(upardim == 1){par_adj[i] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i,])[par_stan])#never happen?
      }else{par_adj[i,] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i,])[par_stan])}
    }#treat 1 dimensional parameter as special due to dimension drop
  }#end else


  #make sure names are the same for sampled and adjusted parms
  row.names(par_adj) <- 1:ndraws
  colnames(par_samps) <- colnames(par_adj)

  if(export_unconst_pars){pars_unc <- list(original = upar_samps, adjusted = upar_adj)}else{pars_unc = NULL}#testing multivariate normality
  if(prior_only){
  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms = par_adj, H = Hhat, J = Jhat, Hprior = H0,
             unconst_parms = pars_unc)
  }else{
    rtn = list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms = par_adj, H = Hhat, J = Jhat, unconst_parms = pars_unc)
  }
  class(rtn) = c("cs_sampling", class(rtn))

  return(rtn)

}#end of cs_sampling_yj
