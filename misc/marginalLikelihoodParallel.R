marginalLikelihoodParallel <- function(mu,
                                     tau2,
                                     Y,
                                     #simulate
                                     sim = 1000
){

    #Set up
    #number of studies
    k = nrow(Y)
    #initialize
    #lik <- 0


    cl <- parallel::makeCluster(4)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("Y", "mu", "tau2"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(stats))

    #Outer loop
    lik <- parallel::parLapply(cl, 1:sim, function(i){

        #Simulate theta and gamma
        #Trouble numerically integrating with so
        #many parameters, so simulate thetas and gammas
        #based on mu and tau2 to optimize and take average lik
        theta <- rnorm(k, mean = mu, sd = sqrt(tau2))
        gamma <- rnorm(k, mean = 0, sd = 10)

        #Calculate odds (expit scale)
        eta_C <- gamma - theta/2
        eta_T <- gamma + theta/2

        #Transform odds to probabilities (p scale)
        p_C <- 1/(1+exp(-eta_C))
        p_T <- 1/(1+exp(-eta_T))

        #Construct priors (on theta and gamma only)
        #Use log scale for easier computation so
        #addition instead of multiplication for
        #likelihood and priors

        pi_theta <- dnorm(theta, mean = mu, sd = sqrt(tau2), log = T)
        pi_gamma <- dnorm(gamma, mean = 0, sd = 10, log = T)

        #Construct likelihood for both control
        #and treatments, because the likelihood
        #depends on the data (do not construct
        #likelihood for the parameters being estimated
        #or integrated out)

        #Use probabilities found above

        lik_C <- dbinom(Y[,"ric"], Y[,"nic"], p_C, log = T)
        lik_T <- dbinom(Y[,"rit"], Y[,"nit"], p_T, log = T)

        #Need to tranform back to exp scale
        return(exp(lik_C+lik_T+pi_theta+pi_gamma))

    })#End outer loop

    lik <- sum(unlist(lik))

    #Find average log likelihood
    log_lik <- log(lik/sim)

    #Return negative log likelihood
    ret <- -sum(log_lik)
    ret

}

marginalLikelihoodParallel(mu = 0.2,
                           tau2 = 0.1,
                           Y = dat,
                           #simulate
                           sim = 100)
