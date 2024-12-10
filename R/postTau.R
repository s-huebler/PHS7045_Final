#' Posterior Tau2
#'
#' @param x A value from a gamma distribution to be transformed
#' @param a_tau Shape
#' @param b_tau Scalar hyperparameter for gamma distribution rate
#' @param k Number of trials
#' @param mu Scalar mean
#' @param thetas Vector of current estimates of length k
#'
#' @return Inverse of a transformed draw from a gamma
#' @export
#'
#' @examples
#' postTau(x = 0.15,
#'         a_tau = 0.001,
#'         k = 7,
#'         b_tau = 0.001,
#'         mu = -0.2,
#'         thetas = c(-0.05, -0.12, 0.20, -0.36, -0.14, -1.03, -0.46))
#'
postTau <- function(x, k, a_tau, b_tau, mu, thetas){
    # #Recall x is originally an inverse of a gamma
     x <- 1/x
    #
    # #Shape does not change so does not require transform
    # #Original shape of gamma distribution should be k/2 + a_tau
    #
    # #Original rate of the gamma distribution should be 1
    # #Transformation given by z*originalrate/newrate
    # new_rate <- 0.5*sum((thetas-mu)^2) + (b_tau)
    # new_gamma <- x*1/new_rate
    #
    # ret <- 1/new_gamma
    # ret

    1/rgamma(1, shape = k/2+a_tau,
             rate = 0.5*sum((thetas-mu)^2)+b_tau)
}
