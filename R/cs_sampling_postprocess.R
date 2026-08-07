#' Adjust draws from an already fitted Stan model
#'
#' `cs_sampling_postprocess` applies the complex-survey adjustment to an existing \code{\link[brms]{brmsfit}} or \code{\link[rstan]{stanfit}} object using \code{\link[csSampling]{cs_sampling}} or, if `yj_transform` is `TRUE`, \code{\link[csSampling]{cs_sampling_yj}}.
#'
#' @param svydes - a \code{\link[survey]{svydesign}} object or a \code{\link[survey]{svrepdesign}} object. This contains cluster ID, strata, and weight information (\code{\link[survey]{svydesign}}) or replicate weight information (\code{\link[survey]{svrepdesign}})
#'
#' @param fit - a \code{\link[brms]{brmsfit}} or \code{\link[rstan]{stanfit}} object
#'
#' @param par_stan - a list of a subset of parameters to output after adjustment. All parameters are adjusted including the derived parameters, so users may want to only compare subsets. The default, NA, will return all parameters.
#'
#' @param data_stan - a list of data inputs used to construct `fit`. This is extracted automatically for a \code{\link[brms]{brmsfit}} object, but is required for a \code{\link[rstan]{stanfit}}.
#'
#' @param rep_design - logical indicating if the svydes object is a \code{\link[survey]{svrepdesign}}. If FALSE, the design will be converted to a \code{\link[survey]{svrepdesign}} using ctrl_rep settings
#'
#' @param ctrl_rep - a list of settings when converting svydes from a \code{\link[survey]{svydesign}} object to a \code{\link[survey]{svrepdesign}} object. replicates - number of replicate weights. type - the type of replicate method to use, the default is mrbbootstrap which sample half of the clusters in each strata to make each replicate (see \code{\link[survey]{as.svrepdesign}}).
#'
#' @param H_estimate - a string indicating the method to use to estimate H. The default "MCMC" is Monte Carlo averaging over posterior draws. Otherwise, a plug-in using the posterior mean.
#'
#' @param matrix_sqrt - a string indicating the method to use to take the "square root" of the R1 and R2 matrices. The default "eigen" uses the eigenvalue decomposition. Otherwise, the Cholesky decomposition is used.
#'
#' @param prior_only - a logical indicating if the stan model has an option for sampling just from the prior distribution. This can be used to further refine the estimates for covariances H and J.
#'
#' @param yj_transform - a logical indicating whether to use the Yeo-Johnson transformation. If `TRUE`, H is estimated from the transformed posterior covariance.
#'
#' @param yj_args - a named list of arguments pass to \code{\link[csSampling]{cs_sampling_yj}}. Including `diag_only`, `subset_matrix`, `export_unconst_pars`, and `yj_range`.
#'
#' @return The output of cs_sampling or cs_sampling_yj.
#'
#' @examples
#'
#' # Survey design information
#' library(survey)
#' data(api)
#' apistrat$wt <- apistrat$pw / mean(apistrat$pw)
#' 
#' dstrat <- svydesign(ids = ~1, strata = ~stype, weights = ~wt, data = apistrat, fpc=~fpc)
#' 
#' # Fit a weighted brms model
#' library(brms)
#'
#' fit <- brm(
#'     bf(api00|weights(wt) ~ ell + meals + mobility, center = FALSE),
#'     data = apistrat,
#'     family = gaussian(),
#'     save_pars = save_pars(all = TRUE)
#' )
#'
#' # Apply the adjustment with Yeo-Johnson transformation 
#' adjusted_yj <- cs_sampling_postprocess(svydes = dstrat, fit = fit)
#'
#' # Apply the adjustment without Yeo-Johnson transformation 
#' adjusted <- cs_sampling_postprocess(
#'     svydes = dstrat, fit = fit, yj_transform = FALSE
#' )
#' 
#' @import rstan 
#' @import brms
#' 
#' @export
#' 
cs_sampling_postprocess <- function(
    svydes, fit, par_stan = NA, data_stan = NULL,
    rep_design = FALSE,
    ctrl_rep = list(replicates = 100, type = "mrbbootstrap"),
    H_estimate = "MCMC",
    matrix_sqrt = "eigen",
    prior_only = FALSE,
    yj_transform = TRUE,
    yj_args = list()
    
) {

    # Get stanfit
    if (inherits(fit, "brmsfit")) {

        if (is.null(data_stan)) {
            data_stan <- brms::standata(fit)
        }

        backend <- fit$backend
        if (is.null(backend)) {
            backend <- "rstan"
        }
        if (!identical(backend, "rstan")) {
            fit <- brms::add_rstan_model(fit) # Convert to rstan style
        }
        stan_fit <- fit$fit

    } else if (inherits(fit, "stanfit")) {
        if (is.null(data_stan)) {
            stop("`data_stan` is required.", call. = FALSE)
        }
        stan_fit <- fit
    } else {
        stop("`fit` must be a `brmsfit` or `stanfit` object.", call. = FALSE)
    }

    if (!inherits(stan_fit, "stanfit")) {
        stop("The fitted object does not contain an RStan `stanfit` object.", call. = FALSE)
    }

    # Get stan mod
    mod_stan <- rstan::get_stanmodel(stan_fit)

    if (yj_transform) {
        yj_call <- list(
            svydes = svydes,
            mod_stan = mod_stan,
            par_stan = par_stan,
            data_stan = data_stan,
            rep_design = rep_design,
            ctrl_rep = ctrl_rep,
            matrix_sqrt = matrix_sqrt,
            prior_only = prior_only,
            stan_fit = stan_fit
        )
        return(do.call(.cs_sampling_yj_process, c(yj_call, yj_args)))
    }

    .cs_sampling_process(
        svydes = svydes,
        mod_stan = mod_stan,
        par_stan = par_stan,
        data_stan = data_stan,
        rep_design = rep_design,
        ctrl_rep = ctrl_rep,
        H_estimate = H_estimate,
        matrix_sqrt = matrix_sqrt,
        prior_only = prior_only,
        stan_fit = stan_fit
    )
}
