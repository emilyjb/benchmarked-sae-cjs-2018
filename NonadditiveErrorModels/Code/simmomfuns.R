
mle.fun <- function(phats){
	mle.out <- glmer(cbind(nis*phats, nis*(1-phats))~ (1|as.factor(1:length(nis))), family = binomial(link = "logit"), nAGQ = 10)
	betahat.mle <- summary(mle.out)$coefficients[,"Estimate"]
	sigma2uhat.mle <- summary(mle.out)$varcor[[1]][[1]]
	c(betahat.mle, sigma2uhat.mle)
}

comp.mom.B <- function(boot.b,B){
	mu <- boot.b[1]
	sigma2u <- boot.b[2]
	phats <- boot.b[-c(1,2)]
	muvec <- rep(mu, m)                                                                                    
	pitv.B <- replicate(B, genpitv.m(muvec, sigma2u))                                                       
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B, phats, nis)
	c(cond.moms[1,], cond.moms[2,] - cond.moms[1,]^2)
}
 	
comp.mom.B.f <- function(boot.b,muvec,sigma2u,B){
	phats <- boot.b[-c(1,2)]
	pitv.B <- replicate(B, genpitv.m(muvec, sigma2u))                                                       
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B, phats, nis)
	c(cond.moms[1,], cond.moms[2,] - cond.moms[1,]^2)
}
 	

genpitv.m  <- function(muvec, sigma2u){
	eta <- muvec + rnorm(length(muvec), 0, 1)*sqrt(sigma2u)
	theta <- exp(eta)/(1 + exp(eta))
	theta
}


condmom.binom.i <- function(k, pitv.B, phats, nis){
	dbinom.k <- sapply(pitv.B[k,], function(x){dbinom(phats[k]*nis[k], size = nis[k], prob = x)})
	p.1.mean <- mean(pitv.B[k,]*dbinom.k)/mean(dbinom.k)
	p.1.mean2 <- mean((pitv.B[k,]^2)*dbinom.k)/mean(dbinom.k)
	c(p.1.mean, p.1.mean2)
}


condvar.binom.i <- function(k, pitv.B, phats, nis){
	dbinom.k <- sapply(pitv.B[k,], function(x){dbinom(phats[k]*nis[k], size = nis[k], prob = x)})
	p.1.mean2 <- mean((pitv.B[k,]^2)*dbinom.k)/mean(dbinom.k)
	p.1.mean2
}

mom.uncond <- function(k, pitv.B, phats, nis){
	dbinom.k <- sapply(pitv.B[k,], function(x){dbinom(phats[k]*nis[k], size = nis[k], prob = x)})
	p.1.mean2 <- mean((pitv.B[k,]^2)*dbinom.k)/mean(dbinom.k)
	p.1.mean2
}


mse.boots <- function(B,B2, muvec, sigma2u){
	thetas.B <- replicate(B, genpitv.m(muvec, sigma2u))
	phats.B <- sapply(1:B, function(b){apply(cbind(nis, thetas.B[,b]), 1, genphat)})
	mle.B <- apply(phats.B, 2, mle.fun)
	boot.B <- rbind(mle.B, phats.B)
	pg1.bstar <- apply(boot.B, 2,comp.mom.B,B=B2)
	pg1.bstar.f <- apply(boot.B, 2,comp.mom.B.f, muvec, sigma2u,B=B2)
	g2.boot <- apply((pg1.bstar[1:length(nis),] - pg1.bstar.f[1:length(nis),])^2,1,mean)
	g1.bias.sub.boot <- apply((pg1.bstar[-c(1:length(nis)),] - pg1.bstar.f[-c(1:length(nis)),]),1,mean)
	g1.bias.div.boot <- apply(pg1.bstar[-c(1:length(nis)),],1,mean)/apply(pg1.bstar.f[-c(1:length(nis)),],1,mean)
	cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot)
}

####  Modify mse.boots.bench to return (1) pg1.bstar, (2) benchmarked boot 
	####  based on linear additive approach 
mse.boots.bench <- function(B,B2, muvec, sigma2u){
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

	g2.bench <- apply((phat.boot.bench.BB - thetas.B)^2, 1, mean)
	
	tstat.abs.BB <- abs((phat.boot.bench.BB - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
	tstat.abs.EB.BB <- abs((pg1.bstar[1:length(nis),] - thetas.B)/sqrt(pg1.bstar[-c(1:length(nis)),]))
	

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
	
	cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2)
}


####  Modify mse.boots.bench to return (1) pg1.bstar, (2) benchmarked boot 
	####  based on linear additive approach 
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
	
	cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2, cstarvec3, g2.aug)
}


####  Modify mse.boots.bench to return (1) pg1.bstar, (2) benchmarked boot 
	####  based on linear additive approach 
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
	
	cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2, cstarvec3, g2.aug)
}



####  Modify mse.boots.bench to return (1) pg1.bstar, (2) benchmarked boot 
	####  based on linear additive approach 
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
	
	cbind(g2.boot, g1.bias.sub.boot, g1.bias.div.boot, g2.bench, cstarvec, cstarvec2, cstarvec3, g2.aug)
}





#####  Calibrate the alpha ..... 
calib.alpha.1 <- function(qstar, tstat.abs.BB){
	tapply(apply(tstat.abs.BB < qstar, 1, mean), nis, mean)[1] - 0.95
}


#####  Calibrate the alpha ..... 
calib.alpha.2 <- function(qstar, tstat.abs.BB){
	tapply(apply(tstat.abs.BB < qstar, 1, mean), nis, mean)[2] - 0.95
}


comp.bootphis <- function(B, pitv.BB){
	ind1s <- (1:B-1)*500 + 1
	ind2s <- ind1s + 499
	phi.boot <- sapply(ind1s, function(i){ apply(pitv.BB[,i:(i+499)]*(1-pitv.BB[,i:(i+499)]),1,mean)/nis + apply(pitv.BB[,i:(i+499)],1,var) })
  	phi.boot
}


comp.laboots <- function(b, phi.inv.BB, phat.EB.BB, phats.BB){
	betahatla.num.b  <- sum(nis*(phats.BB[,b] - phat.EB.BB[,b]))
	betahatla.den.b  <- sum(phi.inv.BB[,b]*nis^2)
	phat.bench.la.b <- phat.EB.BB[,b] + phi.inv*nis*betahatla.num.b/betahatla.den.b
	phat.bench.la.b
}

                                                           
condmean.mult.k <- function(k, Z.gen.M, Z.1, m, nks, K, phatik, Sigma.j, xi.j){
	Pik.B <- genp.ln(sqrtfun(sigma2.u),xi.j,Z.gen.M, Z.1, m,  K)
	Pik.B.k <- Pik.B[c((k-1)*(m-1)+1,k*(m-1)),]
	dmult.k <- apply(Pik.B.k, 2, density.mult, phat = phatik[,k],nk = nks[k])
	p.1.mean <- apply(t(Pik.B.k)*dmult.k, 2, mean)/mean(dmult.k)
	p.condmean <- c(1-sum(p.1.mean), p.1.mean)
}
                                                                       		
condvar.mult.k <- function(k, Z.gen.M, Z.1, m, nks, K, phatik, Sigma.j, xi.j){
	Pik.B <- genp.ln(sqrtfun(Sigma.j),xi.j,Z.gen.M, Z.1, m,  K)
	Pik.B.k <- Pik.B[c((k-1)*(m-1)+1,k*(m-1)),]
	dmult.k <- apply(Pik.B.k, 2, density.mult, phat = phatik[,k],nk = nks[k])
	pcov.B <- cbind(Pik.B.k[1,]*(1-Pik.B.k[1,]), -Pik.B.k[1,]*Pik.B.k[2,], Pik.B.k[2,]*(1-Pik.B.k[2,]))
	p.1.cov <- apply(pcov.B*dmult.k, 2, mean)/mean(dmult.k)
	matrix(c(p.1.cov[1], p.1.cov[2], p.1.cov[2], p.1.cov[3]),2,2)*(nks[k])
}
		


condmean.mult.k <- function(k, Z.gen.M, Z.1, m, nks, K, phatik, Sigma.j, xi.j){
	Pik.B <- genp.ln(sqrtfun(Sigma.j),xi.j,Z.gen.M, Z.1, m,  K)
	Pik.B.k <- Pik.B[c((k-1)*(m-1)+1,k*(m-1)),]
	dmult.k <- apply(Pik.B.k, 2, density.mult, phat = phatik[,k],nk = nks[k])
	p.1.mean <- apply(t(Pik.B.k)*dmult.k, 2, mean)/mean(dmult.k)
	p.condmean <- c(1-sum(p.1.mean), p.1.mean)
}
 