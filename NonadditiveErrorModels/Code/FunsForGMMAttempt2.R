
### Method of moments:

mean.sim.beta <- function(beta, sig2u, Z, n, m, phats){
	eta <- beta + sqrt(sig2u)*Z
	h1 <- n*m*mean(exp(eta)/(1+exp(eta)))
	S <- sum(n*phats)
	dev <- S - h1
	dev
}

mean.sim.sig <- function(lsig2u,beta, Z, n, m, phats){
	eta <- beta + sqrt(exp(lsig2u))*Z
	h2 <- m*n*(n-1)*mean((exp(eta)/(1+exp(eta)))^2)
	S <- sum((n*phats)^2 - (n*phats))
	dev <- S - h2
	dev
}

mean.sim <- function(par, Z, ns, m, phats){
	beta <- par[1]; sig2u <- par[2]
	eta <- beta + sqrt(sig2u)*Z
	h1 <- sum(ns)*mean(exp(eta)/(1+exp(eta)))
	h2 <- sum(ns*(ns-1))*mean((exp(eta)/(1+exp(eta)))^2)
	S1 <- sum(ns*phats)
	S2 <- sum((ns*phats)^2 - (ns*phats))
	#S1 - h1)^2 + (S2 - h2)^2
	c(S1 - h1, S2 - h2)
}

h1fun <- function(par, Z, ns, m, phats){
	beta <- par[1]; sig2u <- par[2]
	eta <- beta + sqrt(sig2u)*Z
	h1 <- sum(ns)*mean(exp(eta)/(1+exp(eta)))
	h1
}
	
h2fun <- function(par, Z, ns, m, phats){
	beta <- par[1]; sig2u <- par[2]
	eta <- beta + sqrt(sig2u)*Z
	h2 <- sum(ns*(ns-1))*mean((exp(eta)/(1+exp(eta)))^2)
	h2
}
	
der.sim <- function(par, Z, n, m, phats){
	beta <- par[1]; sig2u <- par[2]
	alpha <- sqrt(sig2u)*Z
	eta <- beta + alpha
	p <- exp(eta)/(1+exp(eta))
	t11 <- m*n*mean(p*(1-p))
	t12 <- -m*n*mean(p*(alpha^2/4)*(1/sig2u)^3)
	t21 <- m*n*(n-1)*mean(2*p*(1-p))
	t22 <- -m*n*(n-1)*mean(p^2*(alpha^2/4)*(1/sig2u)^3)
	D <- matrix(c(t11, t12, t21, t22), 2, 2, byrow = TRUE)
	D
}

der.sim2 <- function(par, Z, ns, m, phats){
	beta <- par[1]; sigu <- par[2]
	eta <- beta + sigu*Z
	p <- exp(eta)/(1+exp(eta))
	t11 <- sum(ns)*mean(p*(1-p))	
	t12 <- sum(ns)*mean(p*(1-p)*Z)
	t21 <- sum(ns*(ns-1))*mean(2*(p^2)*(1-p))
	t22 <- sum(ns*(ns-1))*mean(2*(p^2)*(1-p)*Z)
	D <- matrix(c(t11, t12, t21, t22), 2, 2, byrow = TRUE)
	D
}

### Function to compute the grid-based initial estimates:
grid.search <- function(seq.mu, seq.sig, fun1, fun2, s1, s2, ns, m, phats){
      Z.temp <- rnorm(1000)
	grid.par <- expand.grid(x=seq.mu, y = seq.sig)
	grid.par$h1 <- apply(grid.par, 1, fun1, Z=Z.temp, ns = ns, m = m, phats = phats)
	grid.par$h2 <- apply(grid.par, 1, fun2, Z=Z.temp, ns = ns, m = m, phats = phats)
	check.0 <- abs(grid.par$h1 - s1) + abs(grid.par$h2 - s2)
	par.init <- grid.par[which.min(check.0),c(1,2)]
	list(par.init, min(check.0))
}

onestep.est <- function(bet.init, sig2u.init, nis, m, phats){
	par.cur <- c(bet.init, sqrt(sig2u.init))
	iter <- 0
	repeat{
		iter <- iter + 1
		Z.temp <- rnorm(1000*m)
		dev.cur <- mean.sim(c(par.cur[1], par.cur[2]^2), Z.temp, nis, m, phats)/m
		D.cur <- der.sim2(par.cur, Z.temp, nis, m, phats)/m
		par.new <- par.cur + solve(D.cur)%*%dev.cur
		maxdiff <- max(abs(par.new - par.cur))
		print(paste(maxdiff))
		if(maxdiff < 10^-3 | iter > 40){break}
		par.cur <- as.vector(par.new)
	}
	c(par.new[1], par.new[2]^2)
}

genpitv.m.f <- function(r, muvec, sigma2u, fixnorm){
        eta <- muvec + fixnorm[,r]*sqrt(sigma2u)
        theta <- exp(eta)/(1 + exp(eta))
        theta
}

eb.aug <- function(delta.param, z.aug, mu.fixed, sigma2uhat.mle, R, fixnorm, phats){ 
	linear.aug <- z.aug*delta.param + mu.fixed
	pitv.B.aug <- sapply(1:R, genpitv.m.f, as.vector(linear.aug), sigma2uhat.mle,fixnorm)
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B.aug, phats, nis)
	(sum(nis*as.vector(cond.moms[1,])) - sum(nis*phats))/sum(nis)
}


eb.aug2 <- function(delta.param, z.aug, mu.fixed, sigma2uhat.mle){ 
	linear.aug <- z.aug*delta.param + mu.fixed
	pitv.B.aug <- sapply(1:R, genpitv.m.f, as.vector(linear.aug), sigma2uhat.mle, fixnorm)
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B.aug, phats, nis)
	cond.moms[1,]
}

reba.b <- function(b, phi.inv.BB, mle.B, fixnorm.list,R, phats.B){
	z.aug.b <- phi.inv.BB[,b]
	mu.fixed.b <- rep(mle.B[1,b], length(nis))
	#fixnorm <- fixnorm.list[[b]]
	root.ebaug.b <- NULL
	try(root.ebaug.b <- uniroot(eb.aug, interval = c(-100, 100), z.aug.b, mu.fixed.b, mle.B[2,b], R, fixnorm.list[[b]], phats.B[,b]))
	if(is.null(root.ebaug.b)){
		rootout <- 0
	}else{
		rootout <- root.ebaug.b$root
	}
	rootout
}

eba2.b <- function(b, roots.B, phi.inv.BB, mle.B, fixnorm.list, R, phats.B){
	z.aug.b <- phi.inv.BB[,b]
	mu.fixed.b <- rep(mle.B[1,b], length(nis))		
	linear.aug.b <- z.aug.b*roots.B[b] + mu.fixed.b
	pitv.B.aug.b <- sapply(1:R, genpitv.m.f, as.vector(linear.aug.b),  mle.B[2,b], fixnorm.list[[b]])
	cond.moms <- sapply(1:length(nis), condmom.binom.i, pitv.B.aug.b, phats.B[,b], nis)
	cond.moms[1,]
}


eb.aug.b <- function(B, phi.inv.BB, mle.B, phats.B, R){
	fixnorm.list <- lapply(as.list(1:B), function(i){matrix(rnorm(D*R), nrow = D)})
	roots.B <- sapply(1:B, reba.b, phi.inv.BB, mle.B, fixnorm.list, R, phats.B)
	phats.aug.B <- sapply(1:B, eba2.b, roots.B, phi.inv.BB, mle.B, fixnorm.list, R, phats.B)
	#pitv.B.aug.b <- sapply(1:R, genpitv.m.f, as.vector(linear.aug.b),  mle.B[2,b], fixnorm.list[[b]])
	phats.aug.B
}


#eb.aug2( -0.07507651, z.aug, mu.fixed, sigma2uhat.mle)





