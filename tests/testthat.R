library(testthat)
library(gmhs2)

test_check("gmhs2")

#Checking that an error should be returned if data arguement does
#not have a column called theta_s
testthat::expect_error(mhGibbs(data.frame("Trial" = c("A", "B"),
                                          "rit" = c(15, 17),
                                          "nit" = c(50, 60),
                                          "ric" = c(13, 13),
                                          "nic" = c(50, 55),
                                          "logOR" = c(0.1, -0.03),
                                          "gamma_i" = c(0.02, 0.02)),
                               hyps, mh_scale = 0.2,
                               sim = 1000, burn = 10,
                               Bayesian = "hierarchical"))
