#' Metropolis-Hastings Single Iteration
#'
#'
#' @param Y Data table with columns rit, nit, ric, nic, theta_i, gamma_i
#' @param mu Estimate for mean
#' @param tau2 Estimate for tau2
#' @param gammas_mat 2xk matrix of current and pregenerated candidate gammas
#' @param thetas_mat 2xk matrix of current and pregenerated candidate thetas
#' @param k The number of trials in Y
#' @param mh_scale MH algorithm variance
#' @param crit Critical value for updating
#' @param hypers Hyper parameters wity priors a_mu, b_mu, a_gamma, b_gamma, a_tau2, b_tau2
#'
#' @return list of new thetas and gammas
#' @export
#'
#' @examples
#' Y <- read.table(text = "
#' trial rit nit ric nic
#' Balcon 14 56 15 58
#' Clausen 18 66 19 64
#' Multicentre 15 100 12 95
#' Barber 10 52 12 47
#' Norris 21 226 24 228
#' Kahler 3 38 6 31
#' Ledwich 2 20 3 20
#' ", header = TRUE)
#' Y$theta_i <- log(Y$rit/Y$nit/(1-Y$rit/Y$nit))-
#' log(Y$ric/Y$nic/(1-Y$ric/Y$nic))
#' Y$gamma_i<- 0.5*(log(Y$rit/Y$nit/(1-Y$rit/Y$nit))+
#' log(Y$ric/Y$nic/(1-Y$ric/Y$nic)))
#' hypers <- list(#Mu priors N(0,100)
#' "a_mu" = 0,
#' "b_mu" = 100,

#' #Gamma priors N(0,100)
#' "a_gamma" = 0,
#' "b_gamma" = 100,

#' #Tau2 priors
#' "a_tau2" = 0.001,
#' "b_tau2" = 0.001)
#'
#' gammas_mat <- matrix(c(-1.0758811, -0.9215264, -1.8342675, -1.2527630,
#' -2.2092769, -1.9419261, -1.9659128,
#' -0.1591379, -1.3580273, -1.4989406,
#' 0.0635526,  0.5562330, -0.5041676, -0.8886074), nrow = 2, ncol = 7,
#' byrow = TRUE)
#' thetas_mat <- matrix(c(-0.04546237, -0.11860574,  0.19933290, -0.36464311,
#' -0.13842138, -1.02961942, -0.46262352,
#' -1.9667471, -0.5626374,  0.3088954, -0.7944453,  1.3856203,  0.3028983,
#'  0.3382107), nrow = 2, ncol = 7,
#' byrow = TRUE)
#' mh(Y, -0.2, 0.07, gammas_mat, thetas_mat, 7, 0.23, -0.3925829, hypers)

mh <- function(Y,      #Data
               mu,     #Current mu estimate
               tau2,   #Current tau2 estimate
               gammas_mat, #Matrix size 2xk
               thetas_mat, #Vector size 2xk
               k,      #Number of observations in meta-analysis
               mh_scale, #Scale of jump
               crit,   #Critical value
               hypers  #Hyperparameters
){

    #Pre allocate memory
    thetas_new <- thetas_mat[1,]
    gammas_new <- gammas_mat[1,]


    #Define candidates (standard normal*jump + previous iter)
    candidates_thetas <- thetas_mat[2,]*mh_scale*sd(thetas_mat[1,])+thetas_mat[1,]
    #print(candidates_thetas)
    candidates_gammas <- gammas_mat[2,]*mh_scale*sd(gammas_mat[1,]) +gammas_mat[1,]



    #Run posterior function on current and proposed data
    # for thetas and gammas and compare to crit value in a single loop


    for(i in 1:k){

        #Thetas
        temp_theta_current <- gmhs2::postTheta(tau2 = tau2,
                                              mu = mu,
                                              gamma = gammas_new[i],
                                              theta = thetas_new[i],
                                              rt = Y[i, "rit"],
                                              rc = Y[i, "ric"]) |>
            log()

       # print(temp_theta_current)

        temp_theta_candidate <- gmhs2::postTheta(tau2 = tau2,
                                                mu = mu,
                                                gamma = candidates_gammas[i],
                                                theta = candidates_thetas[i],
                                                rt = Y[i, "rit"],
                                                rc = Y[i, "ric"]) |>
            log()

       # print(temp_theta_candidate)

        if(is.finite(temp_theta_candidate) && is.finite(temp_theta_current)){
            acceptance_ratio_theta <- exp(temp_theta_candidate - temp_theta_current)
            if(crit < acceptance_ratio_theta){
                thetas_new[i]<- candidates_thetas[i]
            }
        }else{print(paste0("non-finite theta"))}



        #Gammas
        temp_gamma_current <- gmhs2::postGamma(b_gamma = hypers$b_gamma,
                                         tau2 = tau2,
                                         mu = mu,
                                         gamma = gammas_new[i],
                                         theta = thetas_new[i],
                                         rt = Y[i, "rit"],
                                         rc = Y[i, "ric"]) |>
          log()

        temp_gamma_candidate <- gmhs2::postGamma(b_gamma = hypers$b_gamma,
                                           tau2 = tau2,
                                           mu = mu,
                                           gamma = candidates_gammas[i],
                                           theta = candidates_thetas[i],
                                           rt = Y[i, "rit"],
                                           rc = Y[i, "ric"]) |>
          log()


        if(is.finite(temp_gamma_candidate) && is.finite(temp_gamma_current)){
            acceptance_ratio_gamma <- exp(temp_gamma_candidate - temp_gamma_current)
            if(crit < acceptance_ratio_gamma){
                gammas_new[i]<- candidates_gammas[i]
            }
        }else{print(paste0("non-finite gamma"))}


        rm(temp_gamma_candidate,
            temp_gamma_current,
            temp_theta_candidate,
            temp_theta_current)
    } #End for loop

    ret <- list("thetas" = thetas_new,
                "gammas" = gammas_new)

    #print(paste("Thetas=", ret[[1]]))
    #print(paste("Gammas=", ret[[2]]))
    ret
}
