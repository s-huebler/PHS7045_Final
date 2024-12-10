#' Posterior Gamma
#'
#' @param b_gamma Variance for the gamma param
#' @param mu Scalar mu estimate
#' @param tau2 Scalar tau2 estimate
#' @param theta Scalar theta estimate
#' @param gamma Scalar gamma estimate
#' @param rt Scalar successes in the trial group
#' @param rc Scalar successes in the control group
#'
#' @return Updated Gamma
#' @export
#'
#' @examples
#' postGamma(100, 0.2, 0.15, theta = -0.05, gamma = -1.07, rt = 14, rc = 15)
postGamma <- function(b_gamma, mu, tau2, theta, gamma, rt, rc){
    # exp(gamma*(rt-rc))/
    #     ((1+exp(gamma - theta/2))^(rc)*(1+exp(gamma + theta/2))^(rt))*
    #     #1/sqrt(2*pi*b_gamma)*
    #     exp(-1/(2*b_gamma)*(gamma)^2)

    exp(gamma - theta/2)^rc*
        exp(gamma + theta/2)^rt/
        ((1+exp(gamma - theta/2))^(rc)*(1+exp(gamma + theta/2))^(rt))*
        1/sqrt(2*pi*b_gamma)*
        exp(-1/(2*b_gamma)*(gamma)^2)


}
