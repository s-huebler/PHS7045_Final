
#include "RcppArmadillo.h"


// [[Rcpp::depends(RcppArmadillo)]]

/* via the exports attribute we tell Rcpp to make this function
 available from R*/

//' Binomial Coefficient
//'
//' Calculates the binomial coefficient
//'
//' @param n
//' Scalar of sample size
//' @param x
//' Scalar of success count
//' @return
//' Coefficient for binomial pmf
//' @export
// [[Rcpp::export]]
double  binomialCoefficient(int n, int x) {

    //if(x == 0 | x == n){
    //    return 1;
    //} else{
    //    return binomialCoefficient(n-1, x-1)+ binomialCoefficient(n-1, x);
    //}
    if(x > n) return 0;
    double result = 1.0;
    for(int i = 1; i <= x; ++i){
        result *= (n - (x-i))/ static_cast<double>(i);
    }
    return result;
 }

//' Binomial Probability Density
//'
//' Alternative to using dbinom
//'
//' @param n
//' Scalar of sample size
//' @param x
//' Scalar of success
//' @param p
//' Scalar probability of success
//' @return
//' Probability density
//' @export
// [[Rcpp::export]]
double binomialPMF(int n, int x, double p){
     double binomCoeff = binomialCoefficient(n,x);
     return binomCoeff * pow(p,x)*pow(1-p, n-x);
 }



//' Calculates the marginal likelihood
//'
//' Avoiding numeric integration of thetas and gammas
//'
//' @param Y
//' A numeric matrix with columns for ric, nic, rit, nit
//' @param mu
//' The mean
//' @param tau2
//' The variance
//' @param sim
//' The number of simulations
//' @param a_gamma
//' The mean of prior distribution for gammas
//' @param b_gamma
//' The variance of prior distribution for gammas
//' @return
//' Vector of likelihoods
//' @export
// [[Rcpp::export]]
double marginalLikelihood(arma::mat Y, double mu, double tau2, int sim, double a_gamma, double b_gamma) {
     int k = Y.n_rows;
     arma::colvec likelihoods(sim);

     //Loop
     for(int j=0; j < sim; ++j){
        //no pregen
        arma::colvec theta = Rcpp::rnorm(k, mu, sqrt(tau2));
        arma::colvec gamma = Rcpp::rnorm(k, a_gamma, sqrt(b_gamma));

        //Odds
        arma::colvec eta_C = gamma - theta/2;
        arma::colvec eta_T = gamma + theta/2;

        //Probs
        arma::colvec p_C = 1.0/(1.0 + exp(-eta_C));
        arma::colvec p_T = 1.0/(1.0 + exp(-eta_T));

        //Priors
        //arma::colvec pi_theta(k);
        //arma::colvec pi_gamma(k);
        //arma::colvec lik_C(k);
        //arma::colvec lik_T(k);

        //Individual Likelihoods
        arma::colvec trial_likelihoods(k);


        for(int i = 0; i < k; ++i){
            //pi_theta[i] = R::dnorm(theta[i], mu, sqrt(tau2), TRUE);
            double tmp_pi_theta = R::dnorm(theta[i], mu, sqrt(tau2), TRUE);
            //pi_gamma[i] = = R::dnorm(gamma[i], a_gamma, sqrt(b_gamma), TRUE);
            double tmp_pi_gamma = R::dnorm(gamma[i], a_gamma, sqrt(b_gamma), TRUE);

            //lik_C[i]
            double tmp_lik_T = log(binomialPMF(Y(i, 1), Y(i, 0),  p_C[i]));
            //lik_T[i]
            double tmp_lik_C= log(binomialPMF(Y(i, 3), Y(i, 2), p_T[i]));

            trial_likelihoods[i] = tmp_pi_theta + tmp_pi_gamma + tmp_lik_C + tmp_lik_T;
        }

        likelihoods[j] = exp(sum(trial_likelihoods));

     }
     //return(likelihoods);
    //Log of average likelihood
    return log(sum(likelihoods)/sim);
 }


