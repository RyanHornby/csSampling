#' @details 
#' ESTIMATION
#' 
#' \code{\link{cs_sampling}} provides for weighted MCMC estimation of a Bayesian model using survey weights. 
#' It also provides an asymptotic adjustment to the parameter covariance. 
#' It takes in a \code{\link[survey]{svydesign}} object and a \code{\link[rstan]{stan_model}} and inputs for \code{\link[rstan]{sampling}}.
#' \code{\link{cs_sampling_brms}} is a wrapper function that takes inputs in the form of model statements in familiar brms syntax via \code{\link[brms]{brmsformula}}. 
#' The brms helper functions \code{\link[brms]{make_stancode}} and \code{\link[brms]{make_standata}} build Stan models.
#' The wrapper calls call \code{\link{cs_sampling}}.
#' 
#' PLOTTING
#' 
#' \code{\link{cs_sampling}} has an associated plotting method \code{\link{plot.cs_sampling}}, which compares the parameter draws before and after the asymptotic adjustment via
#' \code{\link[GGally]{ggpairs}}. The output can be plotted as a pairs plot or as individual plots.
#' 
#' HELPER FUNCTIONS
#' 
#' \code{\link{cs_sampling}} uses helper functions which manipulate the output from \code{\link[rstan]{sampling}}, 
#' apply a custom estimation function for use with \code{\link[survey]{withReplicates}}, and apply an afine transformation to a 1D array to recenter and scale it.
#' 
#' @examples 
#' library(csSampling)
#' library(rstan)
#' library(brms)
#' library(survey)
#' rstan_options(auto_write = TRUE)
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
#'plot(mod.brms_1)
#subset plot by varnames
#'plot(mod.brms_1, varnames = paste("b", 1:3, sep =""))
#'
#'pp <- plot(mod.brms_1)
#'pp[2,1]
#'
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
#'
#' ####Weighted Logistic Regression for NSDUH####
#' 
#' #import NSDUH data
#' #https://www.datafiles.samhsa.gov/study-dataset/national-survey-drug-use-and-health-2014-nsduh-2014-ds0001-nid16876)#
#' #load(file = "NSDUH_2014.RData")
#'  
#' #subset and clean up the large file#
#' #library(stringr)
#' #change names to all upper case
#' #names(PUF2014_090718) <- str_to_upper(names(PUF2014_090718))
#'
#' #QUESTID2: individual ID
#' #CIGMON: past month smoking
#' #AMDEY2_U: past year depression (adults)
#' #CATAG6: Age groups
#' #ANALWT_C: analysis weights
#' #VESTR: Variance estimation strata
#' #VEREP: Variance estimation PSU (nested within strata)

#' #dat14 <- PUF2014_090718[,c("QUESTID2","CIGMON", "AMDEY2_U", "CATAG6", "ANALWT_C", "VESTR", "VEREP")]
#' #rm(PUF2014_090718);gc(); #clean up memory
#' #Compile Stan model
#' #mod  <- stan_model('wt_logistic.stan')# compile stan code (doi:10.1214/18-BA1143SUPP)
#' 
#' #Skip steps above, dat14 and mod loaded with csSampling
#' 
#' dat14 <- csSampling::dat14
#' mod <- csSampling::mod
#' #subset to adults
#'dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
#'
#' #normalize weights to sum to sample size
#' dat14$WTS <- dat14$ANALWT_C*(length(dat14$ANALWT_C)/sum(dat14$ANALWT_C))
#' #create survey design object#
#' svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)
#'
#' #create list of inputs#
#'X <- model.matrix( ~ AMDEY2_U, data = dat14)
#'y <- dat14$CIGMON
#'k   <- dim(X)[2]
#'n   <- length(y)
#'weights <- dat14$WTS
#'
#'#list of inputs
#'data_stan <- list(y = array(y, dim = n), X = X, k = k, n = n, weights = array(weights, dim = n) )
#'par_stan <- c("beta") #subset of parameters interested in
#'
#'#run Stan and adjustment code
#'set.seed(12345) #set seed to fix output for comparison
#'mod1 <- cs_sampling(svydes = svy14, mod_stan = mod, par_stan = par_stan, 
#'                    data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 40, warmup = 20, thin = 1), sampling_args = list(cores = FALSE))
#'plot(mod1)
#' 

#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
