test_that("Not scaling weights should result in warning", {
  #set up initial data and models (some warnings)
  library(survey)
  library(brms)
  library(rstan)
  dat14 = csSampling::dat14
  dat14$WTS <- dat14$ANALWT_C/mean(dat14$ANALWT_C)
  model_formula <- formula("CIGMON|weights(WTS) ~ AMDEY2_U")
  stancodestr <- make_stancode(brmsformula(model_formula,center = F), data = dat14, family = bernoulli())
  mod = stan_model(model_code = stancodestr)
  #subset to adults
  dat14 <- dat14[as.numeric(dat14$CATAG6) > 1,]
  dat14$WTS <- dat14$ANALWT_C

  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~WTS, data = dat14, nest = TRUE)

  #list of inputs
  data_stan <- make_standata(brmsformula(model_formula,center = F), data = dat14, family = bernoulli())

  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  expect_warning(cs_sampling(svydes = svy14, mod_stan = mod,
                      data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 22, warmup = 20, thin = 1)),
                      "Sum of the weights may not equal the sample size")
})

test_that("Having different weights should result in an error", {
  #normalize weights to sum to sample size
  dat14$WTS <- dat14$ANALWT_C*(length(dat14$ANALWT_C)/sum(dat14$ANALWT_C))

  #create survey design object#
  svy14 <- svydesign(ids = ~VEREP, strata = ~VESTR, weights = ~ANALWT_C, data = dat14, nest = TRUE)

  #update data but not survey design
  model_formula <- formula("CIGMON|weights(WTS) ~ AMDEY2_U")
  data_stan <- make_standata(brmsformula(model_formula,center = F), data = dat14, family = bernoulli())

  #run Stan and adjustment code
  set.seed(12345) #set seed to fix output for comparison
  expect_error(cs_sampling(svydes = svy14, mod_stan = mod,
                           data_stan = data_stan, ctrl_stan = list(chains = 1, iter = 22, warmup = 20, thin = 1)))
})
