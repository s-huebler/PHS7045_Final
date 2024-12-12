dat <-  read.table(text = "
                  trial rit nit ric nic
                 Balcon 14 56 15 58
                Clausen 18 66 19 64
            Multicentre 15 100 12 95
                 Barber 10 52 12 47
                 Norris 21 226 24 228
                 Kahler 3 38 6 31
                Ledwich 2 20 3 20
                  ", header = TRUE)

dat$theta_i <- log(dat$rit/dat$nit/(1-dat$rit/dat$nit))-
    log(dat$ric/dat$nic/(1-dat$ric/dat$nic))

dat$gamma_i <- 0.5*(log(dat$rit/dat$nit/(1-dat$rit/dat$nit))+
                        log(dat$ric/dat$nic/(1-dat$ric/dat$nic)))


