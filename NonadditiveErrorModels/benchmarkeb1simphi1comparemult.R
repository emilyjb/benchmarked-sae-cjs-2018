## Generate data;
rm(list = ls(all = TRUE))
source("Code/SimfunsforChambersmodel1.R")
source("Code/FunsforGMMAttempt2.R")
source("Code/simmomfuns.R")
source("Code/bootfunforcomparemultadj.R")

library("lme4")

mu <- -1.2
sigma2u <- 0.6 
m <- 60
nis <- c(rep(2, m/3), rep(6, m/3), rep(18, m/3))
#nis <- rep(5, m)


phats.store <- c()
thetas.store <- c()
g1hat.store  <- c()
phat.EB.store <- c()
phat.EB.bench.store <- c()
amount2adds <- c()
phat.Bs <- c()
g1hat.tvs <- c()
msehat.boot.adds <- c()
msehat.boot.mults <- c()

mhb.bench.adds <- c()
mhb.bench.mults <- c()

mhb.bench.adds <- c()
mhb.bench.mults <- c()

cmse.benchs <- c()
cmse.ebs <- c()
aug.ebs <- c()

mles <- c()

mhb.aug.adds <- c()
mhb.aug.mults <-c()
cmse.augs <- c()

mhb.bench.mult2s <- c() 
mhb.bench.mult0s <- c()
mhb.aug.mult2s <- c() 
mhb.aug.mult0s <- c()
mhb.bench.mult1s <- c()
mhb.aug.mult1s <- c()



cnt <- 0

R <- 1000
D <- length(nis)


repeat{

	cnt <- cnt + 1

	##### Generates data
	thetas <- replicate(m, genpitv(mu, sigma2u))
	phats <- apply(cbind(nis, thetas), 1, genphat)
	phats.store <- rbind(phats.store, phats)
	thetas.store <- rbind(thetas.store, thetas)
 
  ##### Estimate parameters
	mle.out <- glmer(cbind(nis*phats, nis*(1-phats))~ (1|as.factor(1:length(nis))), family = binomial(link = "logit"), nAGQ = 10)
	betahat.mle <- summary(mle.out)$coefficients[,"Estimate"]
	sigma2uhat.mle <- summary(mle.out)$varcor[[1]][[1]]
	if(sigma2uhat.mle == 0){
		sigma2uhat.mle <- 0.006
	}

	mles <- rbind(mles, c(betahat.mle, sigma2uhat.mle))

	##### Compute EB predictor
	pitv.B <- replicate(500, genpitv.m(rep(betahat.mle, length(nis)), sigma2uhat.mle))
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B, phats, nis)
	phat.EB <- cond.moms[1,]
	g1hat <- cond.moms[2,] - (cond.moms[1,])^2
	phat.EB.store <- rbind(phat.EB.store, phat.EB)
	g1hat.store <- rbind(g1hat.store, g1hat)
	
	##### Compute Phi
	phi.inv <- apply(pitv.B*(1-pitv.B),1,mean)/nis + apply(pitv.B,1,var)
	
	##### Linear additive benchmarking
	betahatla.num <- sum(nis*(phats - phat.EB))
	betahatla.den <- sum(phi.inv*nis^2)
	betahatla <- betahatla.num/betahatla.den
	amount2add <- phi.inv*nis*betahatla

	phat.EB.bench <- phat.EB + amount2add
	phat.EB.bench.store <- rbind(phat.EB.bench.store, phat.EB.bench)


	amount2adds <- rbind(amount2adds, amount2add)
	

	##### Augmented model benchmarking
	z.aug <- phi.inv
	mu.fixed <- rep(betahat.mle, length(nis))
	
	fixnorm <- matrix(rnorm(D*R), nrow = D)
	root.ebaug <- uniroot(eb.aug, interval = c(-100, 100), z.aug, mu.fixed, sigma2uhat.mle, R, fixnorm, phats)
	del <- root.ebaug$root
	aug.eb <- eb.aug2(del, z.aug, mu.fixed, sigma2uhat.mle)
	aug.ebs <- rbind(aug.ebs, aug.eb)

	##### Alternative bootstrap MSE estimators
	msecomp.boot <- mse.boots.bench.withaug(100, 500, rep(betahat.mle,length(nis)), sigma2uhat.mle)
	msehat.boot.add <- g1hat - msecomp.boot[,2] + msecomp.boot[,1]
	msehat.boot.add <- ifelse(msehat.boot.add > 0, msehat.boot.add, g1hat)
	msehat.boot.mult <- g1hat/msecomp.boot[,3] + msecomp.boot[,1]

	msehat.boot.adds <- rbind(msehat.boot.adds, msehat.boot.add)
	msehat.boot.mults <- rbind(msehat.boot.mults, msehat.boot.mult)	

	mhb.bench.add <- g1hat - msecomp.boot[,2] +  msecomp.boot[,1] +  msecomp.boot[,4]
	mhb.bench.add <- ifelse(mhb.bench.add > 0, mhb.bench.add, 0)
	mhb.bench.mult <- g1hat/msecomp.boot[,3] +  msecomp.boot[,1] + msecomp.boot[,4]
	mhb.bench.mult2 <- g1hat^2/msecomp.boot[,9] +  msecomp.boot[,1] + msecomp.boot[,4]
	mhb.bench.mult1 <- g1hat^2/msecomp.boot[,9] +  msecomp.boot[,10]

	mhb.bench.mult0 <- g1hat/msecomp.boot[,3] +  msecomp.boot[,10]

	mhb.bench.mult2s <- rbind(mhb.bench.mult2s, mhb.bench.mult2) 
	mhb.bench.mult1s <- rbind(mhb.bench.mult1s, mhb.bench.mult1) 
	mhb.bench.mult0s <- rbind(mhb.bench.mult0s, mhb.bench.mult0)



	mhb.aug.add <- g1hat - msecomp.boot[,2] +  msecomp.boot[,1] +  msecomp.boot[,8]
	mhb.aug.add <- ifelse(mhb.aug.add > 0, mhb.aug.add, 0)
	mhb.aug.mult <- g1hat/msecomp.boot[,3] +  msecomp.boot[,1] + msecomp.boot[,8]
	mhb.aug.mult2 <- g1hat^2/msecomp.boot[,9] +  msecomp.boot[,1] + msecomp.boot[,8]
	mhb.aug.mult0 <-  g1hat/msecomp.boot[,3]  +  msecomp.boot[,11]
	mhb.aug.mult1 <- g1hat^2/msecomp.boot[,9] +  msecomp.boot[,11]  

	mhb.aug.mult2s <- rbind(mhb.aug.mult2s, mhb.aug.mult2) 
	mhb.aug.mult1s <- rbind(mhb.aug.mult1s, mhb.aug.mult1) 
	mhb.aug.mult0s <- rbind(mhb.aug.mult0s, mhb.aug.mult0)



	mhb.bench.adds <- rbind(mhb.bench.adds, mhb.bench.add)
	mhb.bench.mults <- rbind(mhb.bench.mults, mhb.bench.mult)

	mhb.aug.adds <- rbind(mhb.aug.adds, mhb.aug.add)
	mhb.aug.mults <- rbind(mhb.aug.mults, mhb.aug.mult)
	
	
	cmse.benchs <- rbind(cmse.benchs,(sqrt(g1hat)*msecomp.boot[,5]/1.96)^2)
	cmse.ebs <- rbind(cmse.ebs,(sqrt(g1hat)*msecomp.boot[,6]/1.96)^2)
	cmse.augs <- rbind(cmse.augs,(sqrt(g1hat)*msecomp.boot[,7]/1.96)^2)

	if(cnt%%100 == 0){ save.image("BenchmarkRevisionSimulationUpdate61917.Rdata") }

	print(paste(cnt))
	
	if(cnt == 5000){break}

}

maxcnt <- cnt

mse.emp.eb <- apply((phat.EB.store[1:maxcnt,] - thetas.store[1:maxcnt,])^2,2,mean)
mse.emp.bench <- apply((phat.EB.bench.store[1:maxcnt,] - thetas.store[1:maxcnt,])^2,2,mean)
mse.emp.aug <- apply((aug.ebs[1:maxcnt,] - thetas.store[1:maxcnt,])^2,2,mean)
mse.emp.b <-  apply((phat.Bs[1:maxcnt,] - thetas.store[1:maxcnt,])^2, 2, mean) 


sum(phat.EB.bench.store >1 )
sum(phat.EB.bench.store < 0 )

empmsetab <- cbind( tapply(mse.emp.b, nis, mean),
			  tapply(mse.emp.eb, nis, mean),
		        tapply(mse.emp.bench, nis, mean),
			  tapply(mse.emp.aug, nis, mean))*100

round(empmsetab, 3)


g1.bar <- tapply(apply(g1hat.store , 2, mean), nis, mean)

msebar.add <- tapply(apply(msehat.boot.adds[1:maxcnt,] , 2, mean), nis, mean)
ci.add <- tapply(apply(abs(phat.EB.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(msehat.boot.adds[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msebar.mult <- tapply(apply(msehat.boot.mults[1:maxcnt,] , 2, mean), nis, mean)
ci.mult <- tapply(apply(abs(phat.EB.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(msehat.boot.mults[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msec <- tapply(apply(cmse.ebs , 2, mean), nis, mean)
ci.c <- tapply(apply(abs(phat.EB.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(cmse.ebs[1:maxcnt,]) < 1.96, 2, mean), nis, mean)

msebar.mult.bench <- tapply(apply(mhb.bench.mults[1:maxcnt,], 2, mean), nis, mean)
ci.mult.bench <- tapply(apply(abs(phat.EB.bench.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(mhb.bench.mults[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msebar.add.bench <- tapply(apply(mhb.bench.adds[1:maxcnt,], 2, mean), nis, mean)
ci.add.bench <-  tapply(apply(abs(phat.EB.bench.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(mhb.bench.adds[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msebar.c.bench <- tapply(apply(cmse.benchs[1:maxcnt,], 2, mean), nis, mean)
ci.c.bench <-  tapply(apply(abs(phat.EB.bench.store[1:maxcnt,] - thetas.store[1:maxcnt,])/sqrt(cmse.benchs[1:maxcnt,]) < 1.96, 2, mean), nis, mean)

msebar.add.aug <-  tapply(apply(mhb.aug.adds[1:maxcnt,], 2, mean), nis, mean)
ci.add.aug <- tapply(apply(abs(aug.ebs[1:maxcnt,]- thetas.store[1:maxcnt,])/sqrt(mhb.aug.adds[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msebar.mult.aug <-  tapply(apply(mhb.aug.mults[1:maxcnt,], 2, mean), nis, mean)
ci.mult.aug <- tapply(apply(abs(aug.ebs[1:maxcnt,]- thetas.store[1:maxcnt,])/sqrt(mhb.aug.mults[1:maxcnt,]) < 1.96, 2, mean), nis, mean)
msebar.c.aug <-  tapply(apply( cmse.augs[1:maxcnt,], 2, mean), nis, mean)
ci.c.aug <- tapply(apply(abs(aug.ebs[1:maxcnt,]- thetas.store[1:maxcnt,])/sqrt(cmse.augs[1:maxcnt,]) < 1.96, 2, mean), nis, mean)


outertable <- rbind(cbind(round(cbind(msebar.add, msebar.mult, msec)*100,3), round( cbind(ci.add, ci.mult, ci.c)*100, 1)),
	cbind(round(cbind(msebar.add.bench, msebar.mult.bench, msebar.c.bench)*100,3), round( cbind(ci.add.bench, ci.mult.bench, ci.c.bench)*100, 1)),
	cbind(round(cbind(msebar.add.aug, msebar.mult.aug, msebar.c.aug)*100,3), round( cbind(ci.add.aug, ci.mult.aug, ci.c.aug)*100, 1)))

library("xtable")
xtable(apply(outertable, 2, as.character))



	
