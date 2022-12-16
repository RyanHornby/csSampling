# csSampling

csSampling is an R package with main function cs_sampling.
cs_sampling provides estimation of Bayesian models for data collected from complex survey samples by combining functionality from Stan (via rstan and brms) and the survey package. 
    The user can create a survey weighted model in brms or provide custom weighted model via rstan. Survey design information is provided via svydesign objects from the survey package. 
    cs_sampling estimates the weighted stan model and provides an asymptotic covariance correction for model mis-specification due to using survey sampling weights as plug in values in the likelihood. 
    This is often known as a "design effect" which is the "ratio" between the variance from simple random sample and a complex survey sample.
    See Williams, M. R., and Savitsky, T. D. (2020) Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling. International Statistical Review, https://doi.org/10.1111/insr.12376.
