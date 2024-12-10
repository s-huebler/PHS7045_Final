#' Metropolis-Hastings Gibbs Hybrid Algorithm
#'
#' @import stats
#' @param Y Data matrix
#' @param hypers Vector of hyperparameters
#' @param mh_scale Scalar MH variance scale
#' @param sim Scalar number of simulations
#' @param burn Scalar number of burn in
#' @param Bayesian String either "hierarchical" or "empirical"
#' @param hat_mu Max likelihood for mu
#' @param hat_tau2 Max likelihood for tau2
#'
#' @return Posteriors
#' @export
#'
#' @examples
#'  Y <- read.table(text = "
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
#' hypers <- list(
#' "a_mu" = 0,
#' "b_mu" = 100,
#' "a_gamma" = 0,
#' "b_gamma" = 100,
#' "a_tau2" = 0.001,
#' "b_tau2" = 0.001)
#'
#' mhGibbs(Y, hypers, 0.23, sim = 100, burn = 10, "hierarchical")
#' mhGibbs(Y, hypers, 0.23, sim = 100, burn = 10, "empirical")
mhGibbs<- function(Y,               #Data
                    hypers,  #Hyperparameters
                    mh_scale,              #MH variance scale
                    sim = 100,             #Number of simulations
                    burn = 10,              #Burn in
                    Bayesian = "hierarchical",

                    hat_mu = hypers$a_mu,
                    hat_tau2 = hypers$a_tau2

){

    # Set up
    #number of studies
    k <- nrow(Y)

    # Predraw random variables for transformation
    # Random mus from standard normal
    mu <- stats::rnorm(sim, 0, 1)

    # Random inverse taus from gamma
    tau2 <- 1/stats::rgamma(sim, shape = k/2 + hypers$a_tau2)



    #initial values based on data and hyperparameters
    #thetas initialized by odds ratio in data
    thetas_init <- Y[,"theta_i"]

    #gammas initialized by average logit in data
    gammas_init <- Y[,"gamma_i"]

    #mu initialized by average theta
    mu_init <- mean(thetas_init)

    #tau initialized by variance of thetas
    tau2_init <- sd(thetas_init)^2

    # Matrix of candidate thetas and gammas
    thetas <- matrix(stats::rnorm(k*sim,0, 1), nrow = sim)
    gammas <- matrix(stats::rnorm(k*sim,0, 1), nrow = sim)



    #update these initial values
    mu[1]<- mu_init
    tau2[1]<-tau2_init
    thetas[1,]<- thetas_init
    gammas[1,] <- gammas_init

    #print(thetas)

    # Critical values
    crits <- stats::runif(sim)

    if(Bayesian == "hierarchical"){
        for(j in 2:sim){
            #print(thetas)

            mu[j] <- gmhs2::postMu(x = mu[j-1],
                                  k = k,
                                  sigma = 1/tau2[j-1],
                                  thetas = thetas[j-1,],
                                  b_mu = hypers$b_mu)

            tau2[j] <- gmhs2::postTau(x = tau2[j-1],
                                     b_tau = hypers$b_tau,
                                     mu = mu[j],
                                     thetas = thetas[j-1,])


            temp <- #gmhs2::mh(Y = Y,
                mh(Y=Y,
                             mu = mu[j],
                             tau2 = tau2[j],
                             gammas_mat = gammas[(j-1):j,],
                             thetas_mat = thetas[(j-1):j,],
                             k = k,
                             mh_scale = mh_scale,
                             crit = crits[j],
                             hyper = hypers)


            gammas[j,]<- temp$gammas
            thetas[j,]<- temp$thetas

            rm(temp)


        }
    }# End hierarchical



    if(Bayesian == "empirical"){
        for(j in 2:sim){


            mu[j] <- hat_mu

            tau2[j] <- hat_tau2

            temp <- gmhs2::mh(Y = Y,
                             mu = mu[j],
                             tau2 = tau2[j],
                             gammas_mat = gammas[(j-1):j,],
                             thetas_mat = thetas[(j-1):j,],
                             k = k,
                             mh_scale = mh_scale,
                             crit = crits[j],
                             hyper = hypers)


            gammas[j,]<- temp$gammas
            thetas[j,]<- temp$thetas

            rm(temp)


        }

    }

    # Burn
    ret <-list("mu" = mu[(burn+1):(length(mu))],
               "tau2" = tau2[(burn+1):(length(tau2))],
               "thetas" = thetas[((burn+1):(nrow(thetas))),],
               "gammas" = gammas[((burn+1):(nrow(gammas))),])
    ret
}
