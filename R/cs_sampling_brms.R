#' cs_sampling_brms
#'
#' cs_sampling_brms is a wrapper function that takes inputs in the form of model statements in familiar brms syntax.
#' Then brms helper functions build Stan models and call cs_sampling.
#'
#' @param svydes - a \code{\link[survey]{svydesign}} object or a \code{\link[survey]{svrepdesign}} object. This contains cluster ID, strata, and weight information (\code{\link[survey]{svydesign}}) or replicate weight information (\code{\link[survey]{svrepdesign}})
#' @param brmsmod - \code{\link[brms]{brmsformula}}  object, as input to \code{\link[brms]{make_stancode}}. The \code{\link[brms]{brmsformula}}  must specify a weight variable via weights().
#' @param par_brms - a list of a subset of parameters to output after adjustment. All parameters are adjusted including the derived parameters, so users may want to only compare subsets. The default, NA, will return all parameters.
#' @param data - a data frame, as input to \code{\link[brms]{make_stancode}}
#' @param family - \code{\link[brms]{brmsfamily}} as input to \code{\link[brms]{make_stancode}} specifying distribution and link function
#' @param prior - optional input to \code{\link[brms]{make_stancode}}
#' @param stanvars - optional input to \code{\link[brms]{make_stancode}}
#' @param knots - optional input to \code{\link[brms]{make_stancode}}
#' @param ctrl_stan - a list of control parameters to pass to \code{\link[rstan]{sampling}}. Currently includes the number of chains, iter, warmpup, and thin with defaults.
#' @param rep_design - logical indicating if the svydes object is a \code{\link[survey]{svrepdesign}}. If FALSE, the design will be converted to a \code{\link[survey]{svrepdesign}} using ctrl_rep settings
#' @param ctrl_rep - a list of settings when converting svydes from a \code{\link[survey]{svydesign}} object to a \code{\link[survey]{svrepdesign}} object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate (see \code{\link[survey]{as.svrepdesign}}).
#' @param stancode_args - a list of extra arguments to be passed to \code{\link[brms]{make_stancode}}.
#' @param standata_args - a list of extra arguments to be passed to \code{\link[brms]{make_standata}}.
#' @param H_estimate - a string indicating the method to use to estimate H. The default "MCMC" is Monte Carlo averaging over posterior draws. Otherwise, a plug-in using the posterior mean.
#' @param matrix_sqrt - a string indicating the method to use to take the "square root" of the R1 and R2 matrices. The default "eigen" uses the eigenvalue decomposition. Otherwise, the Cholesky decomposition is used.
#' @param sampling_args - a list of extra arguments that get passed to \code{\link[rstan]{sampling}}.
#' @param sampling_args - a list of extra arguments to be passed to \code{\link[rstan]{sampling}}.
#' @return The output of cs_sampling.
#' @import brms
#' @examples
#'
#' #continuous dependent variable
#' # Survey Design Information
#' library(survey)
#' data(api)
#' apistrat$wt <- apistrat$pw /mean(apistrat$pw)
#'
#' dstrat <-
#' svydesign(id=~1,strata=~stype, weights=~wt, data=apistrat, fpc=~fpc)
#'
#' #Define and Run the Stan Model Via BRMS Wrapper
#' library(brms)
#' set.seed(12345)
#' model_formula <- formula("api00|weights(wt) ~
#'                             ell + meals + mobility")
#' mod.brms <- cs_sampling_brms(svydes = dstrat,
#'                              brmsmod = brmsformula(model_formula, center = FALSE),
#'                              data = apistrat, family = gaussian())
#'
#' #Plot the results
#'
#' plot(mod.brms)
#' plot(mod.brms, varnames = paste("b", 1:4, sep =""))
#'
#' @export
cs_sampling_brms <- function(svydes, brmsmod, data, family, par_brms = NA,prior = NULL, stanvars = NULL, knots = NULL,
                             ctrl_stan = list(chains = 1, iter = 2000, warmup = 1000, thin = 1),
                             rep_design = FALSE, ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
                             stancode_args = list(), standata_args = list(),
                             H_estimate = "MCMC",
                             matrix_sqrt = "eigen",
                             sampling_args = list()) {


  stancode <- do.call(brms::make_stancode, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), stancode_args))
  print("compiling stan model")
  mod_brms  <- rstan::stan_model(model_code = stancode)
  data_brms <- do.call(brms::make_standata, c(list(brmsmod, data = data, family = family, prior = prior, stanvars = stanvars, knots = knots), standata_args))

  return(cs_sampling(svydes = svydes, mod_stan = mod_brms, par_stan = par_brms, data_stan = data_brms,
                     rep_design = rep_design, ctrl_rep = ctrl_rep, ctrl_stan = ctrl_stan,
                     H_estimate = H_estimate, matrix_sqrt = matrix_sqrt, sampling_args))
}

