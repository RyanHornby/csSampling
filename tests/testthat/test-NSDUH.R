

test_that("running with proper inputs should return desired output", {
  
  dat14 = csSampling::dat14
  mod = csSampling::mod
  #subset to adults
  dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
  
  #normalize weights to sum to sample size
  dat14$WTS <- dat14$ANALWT_C*(length(dat14$ANALWT_C)/sum(dat14$ANALWT_C))
  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)
  
  #create list of inputs#
  X <- model.matrix( ~ AMDEY2_U, data = dat14)
  y <- dat14$CIGMON
  k   <- dim(X)[2]
  n   <- length(y)
  weights <- dat14$WTS
  
  #list of inputs
  data_stan <- list(y = array(y, dim = n), X = X, k = k, n = n, weights = array(weights, dim = n) )
  par_stan <- c("beta") #subset of parameters interested in
  
  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  mod1 <- cs_sampling(svydes = svy14, mod_stan = mod, par_stan = par_stan, 
                      data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 40, warmup = 20, thin = 1), sampling_args = list(cores = FALSE))
  
  #mod1$stan_fit@date = ' '
  #attr(mod1$stan_fit@sim$samples[[1]], "elapsed_time") = NULL
  
  #save(mod1, file = "/data/mod1.RData")
  #print("saved")
  
  mod1_expected = csSampling::mod1
  plot_expected = csSampling::plt
  plot_actual   = plot(mod1)
  
  expect_equal(mod1_expected$sampled_parms, mod1$sampled_parms)
  expect_equal(mod1_expected$adjusted_parms, mod1$adjusted_parms)
  expect_equal(plot_actual, plot_expected)
  
})

test_that("Not scaling weights should result in an error", {
  dat14 = csSampling::dat14
  #subset to adults
  dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
  dat14$WTS <- dat14$ANALWT_C
  
  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)
  
  #create list of inputs#
  X <- model.matrix( ~ AMDEY2_U, data = dat14)
  y <- dat14$CIGMON
  k   <- dim(X)[2]
  n   <- length(y)
  weights <- dat14$WTS
  
  #list of inputs
  data_stan <- list(y = array(y, dim = n), X = X, k = k, n = n, weights = array(weights, dim = n) )
  par_stan <- c("beta") #subset of parameters interested in
  
  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  expect_error(cs_sampling(svydes = svy14, mod_stan = mod, par_stan = par_stan, 
                      data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 22, warmup = 20, thin = 1)))
})

test_that("having different weights should result in an error", {
  dat14 = csSampling::dat14
  #subset to adults
  dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
  dat14$WTS <- dat14$ANALWT_C
  
  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)
  #normalize weights to sum to sample size
  dat14$WTS <- dat14$ANALWT_C*(length(dat14$ANALWT_C)/sum(dat14$ANALWT_C))
  
  #create list of inputs#
  X <- model.matrix( ~ AMDEY2_U, data = dat14)
  y <- dat14$CIGMON
  k   <- dim(X)[2]
  n   <- length(y)
  weights <- dat14$WTS
  
  #list of inputs
  data_stan <- list(y = array(y, dim = n), X = X, k = k, n = n, weights = array(weights, dim = n) )
  par_stan <- c("beta") #subset of parameters interested in
  
  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  expect_error(cs_sampling(svydes = svy14, mod_stan = mod, par_stan = par_stan, 
                           data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 22, warmup = 20, thin = 1)))
})

test_that("Having only survey weights should result in an warning", {
  dat14 = csSampling::dat14
  #subset to adults
  dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
  dat14$WTS <- dat14$ANALWT_C
  
  #normalize weights to sum to sample size
  dat14$WTS <- dat14$ANALWT_C*(length(dat14$ANALWT_C)/sum(dat14$ANALWT_C))
  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)
  
  #create list of inputs#
  X <- model.matrix( ~ AMDEY2_U, data = dat14)
  y <- dat14$CIGMON
  k   <- dim(X)[2]
  n   <- length(y)
  weights <- dat14$WTS
  
  #list of inputs
  data_stan <- list(y = array(y, dim = n), X = X, k = k, n = n)
  par_stan <- c("beta") #subset of parameters interested in
  
  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  expect_warning(cs_sampling(svydes = svy14, mod_stan = mod, par_stan = par_stan, 
                           data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 22, warmup = 20, thin = 1)),
                 "No stan data weights, using survey weights instead")
})