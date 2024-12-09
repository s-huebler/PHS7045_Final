#' Posterior Mu
#'
#' @param x A value from a standard normal distribution to be transformed
#' @param k Scalar number of individual trials
#' @param sigma Scalar variance
#' @param thetas Vector of current estimates of length k
#' @param b_mu Scalar hyperparameter for the variance of mu
#'
#' @return A transformed single draw from a standard normal distribution with updated mean and variance
#' @export
#'
#' @examples
#' postMu(x = 0.7,
#'         k = 7,
#'         sigma = 0.19,
#'         thetas = c(-0.05, -0.12, 0.20, -0.36, -0.14, -1.03, -0.46),
#'         b_mu = 100 )
#'
postMu <- function(x, k, sigma, thetas, b_mu){

    new_mean <- sigma*sum(thetas)/(k*sigma + b_mu)
    new_sd <- 1/sqrt(k*sigma + b_mu)

    ret <- new_mean + x * new_sd
    ret
}
