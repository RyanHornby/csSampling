#' load_wt_multi_model
#'
#' @return A string containing a Stan model for weighted multinomial regression
#' @export
load_wt_multi_model <- function(){
 model0 <- paste0("functions{

real wt_multinomial_lpmf(int[,] y, vector lambda, vector weights, int n, int K){
    vector[K] theta;
    real check_term;
	int tmpy[K];
    theta = lambda / sum(lambda);
    check_term  = 0.0;
    for( i in 1:n )
    {
	tmpy = y[i,:];
	check_term    = check_term + weights[i] *  multinomial_lpmf(tmpy | theta);
    }
    return check_term;
  }
} /* end function{} block */

data {
	int<lower=1> n;
	int<lower=0> K;
	int<lower=0, upper = 1> y[n,K];
	vector<lower=0>[n] weights;
	vector<lower=0>[K] alpha;
}

parameters {
  vector<lower=0>[K] lambda;
}
transformed parameters{
  simplex[K] theta = lambda / sum(lambda);
  vector[K] loglam = log(lambda);
}

model {
	//theta ~ dirichlet(alpha);
	lambda    ~ gamma(alpha, 1 );
	target += wt_multinomial_lpmf(y | lambda, weights, n, K);
}"
  )
 return(model0)
}
