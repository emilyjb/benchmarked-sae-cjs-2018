 mse.boots.bench.withaugphi2 <- function(B,B2, muvec, sigma2u){
        thetas.B <- replicate(B, genpitv.m(muvec, sigma2u))
        phats.B <- sapply(1:B, function(b){apply(cbind(nis, thetas.B[,b]), 1, genphat)})
        mle.B <- apply(phats.B, 2, mle.fun)
        boot.B <- rbind(mle.B, phats.B)
        pg1.bstar <- apply(boot.B, 2,comp.mom.B,B=B2)
        #pg1.bstar  <- ifelse(pg1.bstar > 0, pg1.bstar, 0)
        check0 <- sum(pg1.bstar <= 0)
        if(check0 > 0){
                row.means <- apply(pg1.bstar, 1, mean)
                get.rows <- which(apply(pg1.bstar, 1, function(x){ sum(x <= 0)}) >0 )
                get.cols <- which(apply(pg1.bstar, 2, function(x){ sum(x <= 0)}) >0 )
                pg1.bstar[get.rows,get.cols] <- row.means[get.rows]
        }
        pg1.bstar.f <- apply(boot.B, 2,comp.mom.B.f, muvec, sigma2u,B=B2)
        g2.boot <- apply((pg1.bstar[1:length(nis),] - pg1.bstar.f[1:length(nis),])^2,1,mean)
        g1.bias.sub.boot <- apply((pg1.bstar[-c(1:length(nis)),] - pg1.bstar.f[-c(1:length(nis)),]),1,mean)
        g1.bias.div.boot <- apply(pg1.bstar[-c(1:length(nis)),],1,mean)/apply(pg1.bstar.f[-c(1:length(nis)),],1,mean)

        pitv.BB <- replicate(500*B, genpitv.m(muvec, sigma2u))
        phi.inv.BB <- pg1.bstar[-c(1:length(nis)),]
        phat.boot.bench.BB <- sapply(1:B, comp.laboots, phi.inv.BB, pg1.bstar[1:length(nis),], phats.B)
        phats.aug.BB <- eb.aug.b(B, phi.inv.BB, mle.B, phats.B, R)

        g2.bench <- apply((phat.boot.bench.BB - pg1.bstar[1:length(nis),] )^2, 1, mean)
        g2.aug <- apply((phats.aug.BB - pg1.bstar[1:length(nis),] )^2, 1, mean) 
        g2.bench2 <- apply((phat.boot.bench.BB - pg1.bstar.f[1:length(nis),] )^2, 1, mean)
        g2.aug2 <- apply((phats.aug.BB - pg1.bstar.f[1:length(nis),] )^2, 1, mean)

        tstat.abs.BB <- abs((phat.boot.bench.BB - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
        tstat.abs.EB.BB <- abs((pg1.bstar[1:length(nis),] - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
        tstat.abs.aug <- abs((phats.aug.BB - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))

        cstar.1 <- cstar.2 <- NA        
        try(cstar.1 <- uniroot(calib.alpha.1, interval = c(0 , 10), tstat.abs.BB)$root)
        if(is.na(cstar.1)){
                cstar.1 <- 10
        }
        try(cstar.2 <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.BB)$root)
        if(is.na(cstar.2)){
                cstar.2 <- 10
        }

        cstar.1.EB <- cstar.2.EB <- NA
        try(cstar.1.EB <- uniroot(calib.alpha.1, interval = c(0, 10), tstat.abs.EB.BB)$root)
        if(is.na(cstar.1.EB)){
                cstar.1.EB <- 10
        }
        try(cstar.2.EB <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.EB.BB)$root)
        if(is.na(cstar.2.EB)){
                cstar.2.EB <- 10
        }

        cstar.1.aug <- cstar.2.aug <- NA
        try(cstar.1.aug <- uniroot(calib.alpha.1, interval = c(0, 10), tstat.abs.aug)$root)
        if(is.na(cstar.1.aug)){
                cstar.1.aug <- 10
        }
        try(cstar.2.aug <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.aug)$root)
        if(is.na(cstar.2.aug)){
                cstar.2.aug <- 10
        }

	  g1.bias.div.boot.2 <- apply(pg1.bstar[-c(1:length(nis)),],1,mean)
	

      # phi.inv.BB <- apply(pitv.BB*(1-pitv.BB),1,mean)/nis + apply(pitv.BB,1,var)
        #betahatla.num <- sum(nis*(phats - phat.EB))
        #betahatla.den <- sum(phi.inv*nis^2)

        cstarvec <- rep(0,length(nis))
        cstarvec[nis == 2] <- cstar.1
        cstarvec[nis == 6] <- cstar.2
        cstarvec[nis == 18] <- cstar.2

        cstarvec2 <- rep(0,length(nis))
        cstarvec2[nis == 2] <- cstar.1.EB
        cstarvec2[nis == 6] <- cstar.2.EB
        cstarvec2[nis == 18] <- cstar.2.EB

        cstarvec3 <- rep(0,length(nis))
        cstarvec3[nis == 2] <- cstar.1.aug
        cstarvec3[nis == 6] <- cstar.2.aug
        cstarvec3[nis == 18] <- cstar.2.aug
        
        cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2, cstarvec3, g2.aug, g1.bias.div.boot.2, g2.bench2, g2.aug2)
}

mse.boots.bench.withaug <- function(B,B2, muvec, sigma2u){
        thetas.B <- replicate(B, genpitv.m(muvec, sigma2u))
        phats.B <- sapply(1:B, function(b){apply(cbind(nis, thetas.B[,b]), 1, genphat)})
        mle.B <- apply(phats.B, 2, mle.fun)
        boot.B <- rbind(mle.B, phats.B)
        pg1.bstar <- apply(boot.B, 2,comp.mom.B,B=B2)
        #pg1.bstar  <- ifelse(pg1.bstar > 0, pg1.bstar, 0)
        check0 <- sum(pg1.bstar <= 0)
        if(check0 > 0){
                row.means <- apply(pg1.bstar, 1, mean)
                get.rows <- which(apply(pg1.bstar, 1, function(x){ sum(x <= 0)}) >0 )
                get.cols <- which(apply(pg1.bstar, 2, function(x){ sum(x <= 0)}) >0 )
                pg1.bstar[get.rows,get.cols] <- row.means[get.rows]
        }
        pg1.bstar.f <- apply(boot.B, 2,comp.mom.B.f, muvec, sigma2u,B=B2)
        g2.boot <- apply((pg1.bstar[1:length(nis),] - pg1.bstar.f[1:length(nis),])^2,1,mean)
        g1.bias.sub.boot <- apply((pg1.bstar[-c(1:length(nis)),] - pg1.bstar.f[-c(1:length(nis)),]),1,mean)
        g1.bias.div.boot <- apply(pg1.bstar[-c(1:length(nis)),],1,mean)/apply(pg1.bstar.f[-c(1:length(nis)),],1,mean)

        pitv.BB <- replicate(500*B, genpitv.m(muvec, sigma2u))
        phi.inv.BB <- comp.bootphis(B, pitv.BB)
        phat.boot.bench.BB <- sapply(1:B, comp.laboots, phi.inv.BB, pg1.bstar[1:length(nis),], phats.B)
        phats.aug.BB <- eb.aug.b(B, phi.inv.BB, mle.B, phats.B, R)

        g2.bench <- apply((phat.boot.bench.BB - pg1.bstar[1:length(nis),] )^2, 1, mean)
        g2.aug <- apply((phats.aug.BB - pg1.bstar[1:length(nis),] )^2, 1, mean) 
        g2.bench2 <- apply((phat.boot.bench.BB - pg1.bstar.f[1:length(nis),] )^2, 1, mean)
        g2.aug2 <- apply((phats.aug.BB - pg1.bstar.f[1:length(nis),] )^2, 1, mean)

        tstat.abs.BB <- abs((phat.boot.bench.BB - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
        tstat.abs.EB.BB <- abs((pg1.bstar[1:length(nis),] - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
        tstat.abs.aug <- abs((phats.aug.BB - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))

        cstar.1 <- cstar.2 <- NA        
        try(cstar.1 <- uniroot(calib.alpha.1, interval = c(0 , 10), tstat.abs.BB)$root)
        if(is.na(cstar.1)){
                cstar.1 <- 10
        }
        try(cstar.2 <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.BB)$root)
        if(is.na(cstar.2)){
                cstar.2 <- 10
        }

        cstar.1.EB <- cstar.2.EB <- NA
        try(cstar.1.EB <- uniroot(calib.alpha.1, interval = c(0, 10), tstat.abs.EB.BB)$root)
        if(is.na(cstar.1.EB)){
                cstar.1.EB <- 10
        }
        try(cstar.2.EB <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.EB.BB)$root)
        if(is.na(cstar.2.EB)){
                cstar.2.EB <- 10
        }

        cstar.1.aug <- cstar.2.aug <- NA
        try(cstar.1.aug <- uniroot(calib.alpha.1, interval = c(0, 10), tstat.abs.aug)$root)
        if(is.na(cstar.1.aug)){
                cstar.1.aug <- 10
        }
        try(cstar.2.aug <- uniroot(calib.alpha.2, interval = c(0 , 10), tstat.abs.aug)$root)
        if(is.na(cstar.2.aug)){
                cstar.2.aug <- 10
        }

	  g1.bias.div.boot.2 <- apply(pg1.bstar[-c(1:length(nis)),],1,mean)
	

    #  phi.inv.BB <- apply(pitv.BB*(1-pitv.BB),1,mean)/nis + apply(pitv.BB,1,var)
        #betahatla.num <- sum(nis*(phats - phat.EB))
        #betahatla.den <- sum(phi.inv*nis^2)

        cstarvec <- rep(0,length(nis))
        cstarvec[nis == 2] <- cstar.1
        cstarvec[nis == 6] <- cstar.2
        cstarvec[nis == 18] <- cstar.2

        cstarvec2 <- rep(0,length(nis))
        cstarvec2[nis == 2] <- cstar.1.EB
        cstarvec2[nis == 6] <- cstar.2.EB
        cstarvec2[nis == 18] <- cstar.2.EB

        cstarvec3 <- rep(0,length(nis))
        cstarvec3[nis == 2] <- cstar.1.aug
        cstarvec3[nis == 6] <- cstar.2.aug
        cstarvec3[nis == 18] <- cstar.2.aug
        
        cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2, cstarvec3, g2.aug, g1.bias.div.boot.2, g2.bench2, g2.aug2)
}





