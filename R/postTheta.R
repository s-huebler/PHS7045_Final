#' Posterior Theta
#'
#' @param tau2 Current tau2 estimate
#' @param mu Current mu estimate
#' @param gamma Scalar gamma
#' @param theta Scalar gamma estimate
#' @param rt Scalar successes in the trial group
#' @param rc Scalar successes in the control group
#'
#' @return Updated Gamma
#' @export
#'
#' @examples
#' postTheta(0.07, -0.2, gamma = -1.34, theta = -0.05, rt = 14, rc = 15)
#'
#'
postTheta <- function(tau2, mu, gamma, theta, rt, rc){
    exp(theta/2*(rt-rc))/
        ((1+exp(gamma - theta/2))^(rc)*(1+exp(gamma + theta/2))^(rt))*
        1/sqrt(2*pi*tau2)*
        exp(-1/(2*tau2)*(theta-mu)^2)

    # exp(gamma - theta/2)^rc*
    #     exp(gamma + theta/2)^rt/
    #     ((1+exp(gamma - theta/2))^(rc)*(1+exp(gamma + theta/2))^(rt))*
    #     1/sqrt(2*pi*tau2)*
    #     exp(-1/(2*tau2)*(theta-mu)^2)
}
