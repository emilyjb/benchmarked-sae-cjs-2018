
genpitv <- function(mu, sigma2u){
	eta <- mu + rnorm(1)*sqrt(sigma2u)
	theta <- exp(eta)/(1 + exp(eta))
	theta
}

genphat <- function(nipi){
	ni <- nipi[1]; pi <- nipi[2]
	rbinom(1, size = ni, prob = pi)/ni
}

iter.beta <- function(sig2uhat.cur, u.cur, sig2ehats, theta.cur){
	derg <- 1/(theta.cur*(1-theta.cur))
	sig2ehats <- (theta.cur*(1-theta.cur))/nis
	V.cur <- sig2uhat.cur + derg*sig2ehats*derg
	muhat.update <- sum(u.cur/V.cur)/sum(1/V.cur)
	muhat.update
}

iter.b <- function(sig2uhat.cur, u.cur, sig2ehats, muhat.cur, theta.cur){
	derg <- 1/(theta.cur*(1-theta.cur))
	sig2ehats <- (theta.cur*(1-theta.cur))/nis
	V.cur <- sig2uhat.cur + derg*sig2ehats*derg
	b.new <- sig2uhat.cur/V.cur*(u.cur - muhat.cur)
	b.new
}

iter.sigma <- function(b.cur, sig2u.cur, sig2ehats, theta.cur){
	derg <- 1/(theta.cur*(1-theta.cur))
	sig2ehats <- (theta.cur*(1-theta.cur))/nis
	V.cur <- sig2uhat.cur + derg*sig2ehats*derg
	term1 <- mean(b.cur^2)
	term2 <- mean(1/(1/V.cur + 1/sig2uhat.cur))
	term1 + term2
}

iter.all.onestep <- function(sig2u, u, sig2e, muhat, b, theta){
	muhat.new <-  iter.beta(sig2u, u, sig2e, theta)
	b.new <- iter.b(sig2u, u, sig2e, muhat, theta)
	sig2uhat.new <- iter.sigma(b, sig2u, sig2e, theta)
	list(muhat.new, b.new, sig2uhat.new)
}


iter.all <- function(sig2u, u, sig2e, muhat, b, theta, phats){
	params.cur <- iter.all.onestep(sig2u, u, sig2e, muhat, b, theta)
	iter <- 0
	repeat{
		iter <- iter + 1
		muhat.cur <- params.cur[[1]]
		theta.cur <- exp(muhat.cur + params.cur[[2]])/(1+exp(muhat.cur + params.cur[[2]]))
		u.cur <- params.cur[[1]] + params.cur[[2]] +  1/(theta.cur*(1-theta.cur))*(phats - theta.cur)
		params.new <- iter.all.onestep(params.cur[[3]], u.cur, sig2ehats, params.cur[[1]], params.cur[[2]], theta.cur)
		diff.musig <- c(params.cur[[1]] - params.new[[1]], params.cur[[2]] - params.new[[2]])
		max.diff <- max(abs(diff.musig))
		params.cur <- params.new
		print(paste(iter))
		print(paste(max.diff))
		if(max.diff < 10^(-3) | iter == 100){break}
	}
	muhat.cur <- params.cur[[1]]
	theta.cur <- exp(muhat.cur + params.cur[[2]])/(1+exp(muhat.cur + params.cur[[2]]))
	derg <- 1/(theta.cur*(1-theta.cur))
	sig2ehats <- (theta.cur*(1-theta.cur))/nis
	V.cur <- params.cur[[3]] + derg*sig2ehats*derg
	params.cur[[4]] <- V.cur
	params.cur
	
}


	

	





