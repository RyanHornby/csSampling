#' cs_sampling
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
#' @param H_estimate - a string indicating the method to use to estimate H. The default "MCMC" is Monte Carlo averaging over posterior draws. Otherwise, a plug-in using the posterior mean.
#'
#' @param matrix_sqrt - a string indicating the method to use to take the "square root" of the R1 and R2 matrices. The default "eigen" uses the eigenvalue decomposition. Otherwise, the Cholesky decomposition is used.
#'
#' @param sampling_args - a list of extra arguments that get passed to \code{\link[rstan]{sampling}}.
#'
#' @import rstan
#' @import survey
#' @import plyr
#' @import pkgcond
#'
#'
#' @examples
#' # multinomial dependent variable
#'
#' library(survey)
#' data(api)
#' #make sure weights sum to n
#' apiclus1$newpw <- apiclus1$pw/mean(apiclus1$pw)
#' dclus1<-svydesign(id=~dnum, weights=~newpw, data=apiclus1, fpc=~fpc)
#' M1 <- svymean(~stype, dclus1)
#'
#' #Use default replicate design
#' rclus1<-as.svrepdesign(dclus1)
#' M2 <- svymean(~stype, rclus1)
#'
#' #Construct a Stan model
#' library(rstan)
#' model_string <- csSampling::load_wt_multi_model()
#' mod_dm <- stan_model(model_code = model_string)
#'
#' #Set the Data for Stan
#' y <- as.factor(rclus1$variables$stype)
#' yM <- model.matrix(~y -1)
#' n <- dim(yM)[1]
#' K <- dim(yM)[2]
#' alpha<-rep(1,K)
#' weights <- rclus1$pweights
#' data_stan<-list("y"=yM,"alpha"=alpha,"K"=K, "n" = n, "weights" = weights)
#' ctrl_stan<-list("chains"=1,"iter"=2000,"warmup"=1000,"thin"=1)
#'
#' #Run the Stan mode
#' set.seed(12345)
#' mod1 <- cs_sampling(svydes = rclus1, mod_stan = mod_dm,
#' data_stan = data_stan, ctrl_stan = ctrl_stan, rep_design = TRUE)
#'
#' #Plot the results
#'
#' plot(mod1)
#' plot(mod1, varnames = paste("theta",1:3, sep = ""))
#'
#' #Compare Summary Statistics
#'
#' #hybrid variance adjusted
#' cbind(colMeans(mod1$adjusted_parms[,4:6]),
#' sqrt(diag(cov(mod1$adjusted_parms[,4:6]))))
#' #Taylor Linearization
#' M1
#' #Replication
#' M2
#'
#' #Bernoulli-distributed dependent variable
#'
#' dat14 <- csSampling::dat14
#' #subset to adults
#' dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
#' #normalize weights to sum to sample size
#' dat14$WTS <- dat14$ANALWT_C/mean(dat14$ANALWT_C)
#' svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR,
#'                    weights = ~WTS, data = dat14, nest = TRUE)
#'
#' #Construct a Stan model
#' library(brms)
#' model_formula <- formula("CIGMON|weights(WTS) ~ AMDEY2_U")
#' stancode <- make_stancode(brmsformula(model_formula,center = FALSE),
#'                           data = dat14, family = bernoulli(), save_model = "brms_wt_log.stan")
#' mod_brms  <- stan_model('brms_wt_log.stan')
#'
#' #Prepare data for Stan modelling
#'
#' data_brms <- make_standata(brmsformula(model_formula,center = FALSE),
#'                           data = dat14, family = bernoulli())
#'
#' #Run the Stan model
#'
#' set.seed(12345) #set seed to fix output for comparison
#' mod.brms <- cs_sampling(svydes = svy14, mod_stan = mod_brms, data_stan = data_brms)
#'
#'
#' #Plot the results
#'
#' plot(mod.brms, varnames = paste("b", 1:2, sep =""))
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
cs_sampling <- function(svydes, mod_stan, par_stan = NA, data_stan,
                        ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                        rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                        H_estimate = "MCMC",
                        matrix_sqrt = "eigen",
                        sampling_args = list()){

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


  print("(1) stan fitting (1)")
  out_stan  <- do.call(rstan::sampling, c(list(object = mod_stan, data = data_stan,
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

  #number of MCMC draws
  ndraws <- dim(par_samps)[1]

  print("(2) Transforming Parameters (2)")
  #convert to list type input > convert to unconstrained parameterization > back to matrix/array
  ##This is the bottleneck###

  #preallocate upar_samps ahead of time instead of using rbind
  tmplist <- list_2D_row_subset(par_samps_list, 1)
  upar_samps_init <- unconstrain_pars(out_stan, tmplist)
  upar_samps <- matrix(data = NA, nrow = ndraws, ncol = length(upar_samps_init))

  if(H_estimate == "MCMC"){ #Average Hessian across MCMC draws
    for (i in 1:ndraws) {
      if(i %% 500 == 0){print(paste0("Converting draw ", i))} #status message for user
      tmplist <- list_2D_row_subset(par_samps_list, i)
      upar_samps[i,] <- rstan::unconstrain_pars(out_stan, tmplist)
      # just need to convert rvar type to numeric within each variable (list entry)
      if (i == 1) {
      #add H estimates
        Hmcmc <- -1 * stats::optimHess(upar_samps[i,], gr = function(x) {grad_log_prob(out_stan, x)})/ndraws
      }
      else {
        Hmcmc <- Hmcmc  - 1 * stats::optimHess(upar_samps[i,], gr = function(x) {grad_log_prob(out_stan, x)})/ndraws
      }
    }
  }else{
    for(i in 1:ndraws){#just need the length here
      if(i %% 500 == 0){print(paste0("Converting draw ", i))} #status message for user
      tmplist <- list_2D_row_subset(par_samps_list, i)
      upar_samps[i,] <- rstan::unconstrain_pars(out_stan, tmplist)
    }
  }

  row.names(upar_samps) <- 1:dim(par_samps)[1]

  upar_hat <- colMeans(upar_samps)

  #Estimate Hessian
  if(H_estimate == "MCMC"){
    Hhat <- Hmcmc
  }else{#use posterior mean plug-in
    Hhat  <- -1*stats::optimHess(upar_hat, gr = function(x){rstan::grad_log_prob(out_stan, x)})
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
  Jhat <- stats::vcov(rep_tmp)

  print("(4) Estimating Adjustment (4)")
  #compute adjustment
  #use pivot for numerical stability - close to positive semi-definite if some parameters are highly correlated
  #(Q <- chol(m, pivot = TRUE))
  ## we can use this by
  #pivot <- attr(Q, "pivot")
  #Q[, order(pivot)]
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
  upar_adj <- plyr::aaply(upar_samps, 1, DEadj, par_hat = upar_hat, R2R1 = R2iR1, .drop = TRUE)

  #back transform to constrained parameter space
  #special cases (4) needed for combinations of constrained and unconstrained par have 1 dimension
  #not likely/possible for unconstrained dim > constrained so one case might never be used

  print("(6) Backtransforming Parameters (6)")
  if(is.null(dim(upar_adj))){upardim <- length(upar_adj)
    upardim <- dim(upar_adj)[1] #treat 1 dimensional parameter as special due to dimension drop
    par_adj_tmp <- unlist(rstan::constrain_pars(out_stan, upar_adj[1])[par_stan])
    par_adj <- matrix(data = NA, nrow = upardim, ncol = length(par_adj_tmp))
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
    for (i in 1:upardim) {
      if(i %% 500 == 0){print(paste0("Back Converting draw ", i))} #status message for user
      if(upardim == 1){par_adj[i] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i,])[par_stan])#never happen?
      }else{par_adj[i,] <- unlist(rstan::constrain_pars(out_stan, upar_adj[i,])[par_stan])}
    }#treat 1 dimensional parameter as special due to dimension drop
  }#end else


  #make sure names are the same for sampled and adjusted parms
  row.names(par_adj) <- 1:ndraws
  colnames(par_samps) <- colnames(par_adj)

  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))

  return(rtn)

}#end of cs_sampling
