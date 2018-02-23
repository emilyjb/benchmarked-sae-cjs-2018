

genudir <- function(provprop, kappa, phi){
  
  alpha0 <- (1-phi)^2/kappa - 1
  alphas <- provprop*alpha0
  gammas <- sapply(alphas, rgamma, n=1, scale=1)
  dirs <- gammas/sum(gammas)
  
  dirs
  
}

gendatafun <- function(i, censustab, Tdotk, nks, ncs, deltasmat, alphasvec, taus, fakeck, K){
  
  genocslist.temp <- vector("list", 4)
  
  pstar <- apply(censustab, 2, genudir, kappa=kappa, phi=0)
  Mstar <- t(t(pstar)*Tdotk)
  uik <- pstar - censustab
  
  ##Generate sample design proportions
  MhatikMtildedes <- genmixeddirsnozerosspNEWb(pstar, Tdotk, nks, ncs,deltasmat,alphasvec,taus,fakeck, K)
    
  genocslist.temp$phatik <- MhatikMtildedes[[1]]
  genocslist.temp$pstar <-pstar
  genocslist.temp$Pik <- censustab
  genocslist.temp$Sighateedir <- MhatikMtildedes[[2]]
  genocslist.temp$nks <- nks
  genocslist.temp$Tdotk <- Tdotk
  genocslist.temp
  
}


genmixeddirsnozerosspNEWb <- function(pstar, Tdotk, nks, ncs,deltasmat,alphasvec,taus,fakeck, K){

m <- dim(pstar)[1]
out1 <- lapply(1:K, gendirfrompssp, pstar,deltasmat, alphasvec, taus, ncs, nks)
out1[[1]] <- out1[[1]]%*%kronecker(diag(nks[1]/2), rep(1,2))
out1[[8]] <- out1[[8]]%*%kronecker(diag(nks[8]/2), rep(1,2))
ncsb <- ncs
ncsb[1] <- 2; ncsb[8] <- 2
phatlist <- lapply(as.list(1:K), getmeanfromclustsp, out1,ncsb)
sigelist <- lapply(as.list(1:K), getvarfromclustsp, out1,ncsb)
Sighatee <- as.matrix(bdiag(sigelist))
phatik <- matrix(unlist(phatlist),m,K)
list(phatik, Sighatee)

}

gendirfrompssp <- function(k, Pik,deltasmat, alphasvec, taus, ncs, nks){
  rksk <- nks[k]/ncs[k]
  replicate(rksk,genmixeddirsfromps1sp(k, Pik,deltasmat, alphasvec, taus, ncs), simplify = TRUE)
}

genmixeddirsfromps1sp <- function(k, Pik,deltasmat, alphasvec, taus, ncs){
  
  betap <- getbetasfrompsp(k, Pik, deltasmat, alphasvec)
  psiml <- getpsfrombetas(betap)
  if(!is.na(taus[k])){
    pclust <- rmultinom(n=1,size=(taus[k] + 1), prob=psiml)/(taus[k]+1)
  }else{
    pclust <- Pik[,k]
  }
  #pclust <- genudir(psiml, 1/(1+taus[k]), 0)
  pmultclus <- rmultinom(pclust, n=1, size=ncs[k])
  pmultclus
  
}


getbetasfrompsp <- function(k, ptemp,deltas,alphas){
  p1.1 <- ptemp[1:(m-1),k]/deltas[,k]
  p1m <- 1-sum(p1.1)
  b1 <- (ptemp[m,k]/p1m)*alphas[k]
  p1 <- c(p1.1, p1m)
  p2 <- (ptemp[,k] - b1*p1)/(1-b1)
  #p2*(1-b1) + p1*b1
  list(b1, cbind(p1, p2))
  
}



getpsfrombetas <- function(betaps){
  if(betaps[[1]]<1){
    dind1 <- genbinom(c(1, betaps[[1]]))
    dind <- c(1,0)*dind1 + c(0,1)*(1-dind1)
    psimnew <- betaps[[2]]%*%dind
    psimnew
  }else{
    betaps[[2]][,1]
  }
  
}


genbinom <- function(dfnp){
  n <- dfnp[1] ; p <- dfnp[2]
  us <- runif(n)
  b <- sum(us<=p)
  b
}


getmeanfromclustsp <- function(k,multlist,ncs){
  mult <- multlist[[k]]; nc <- ncs[k]
  phat <- apply(mult, 1, mean)/nc
  phat
}

getvarfromclustsp <- function(k,multlist,ncs){
  mult <- multlist[[k]]; nc <- ncs[k]
  r <- dim(mult)[2]
  vhat <- var(t(mult))/nc/r/nc
  vhat
}



genmixeddirsnozerosspNEWbNRI <- function(pstar, Tdotk, nks, ncs,deltasmat,alphasvec,taus,fakeck, K){

m <- dim(pstar)[1]
out1 <- lapply(1:K, gendirfrompssp, pstar,deltasmat, alphasvec, taus, ncs, nks)
ncsb <- ncs
#ncsb[1] <- 2; ncsb[8] <- 2
phatlist <- lapply(as.list(1:K), getmeanfromclustsp, out1,ncsb)
sigelist <- lapply(as.list(1:K), getvarfromclustsp, out1,ncsb)
Sighatee <- as.matrix(bdiag(sigelist))
phatik <- matrix(unlist(phatlist),m,K)
list(phatik, Sighatee)

}



gendatafunNRI <- function(i, censustab, Tdotk, nks, ncs, deltasmat, alphasvec, taus, fakeck, K){
  
  genocslist.temp <- vector("list", 4)
  
  pstar <- apply(censustab, 2, genudir, kappa=kappa, phi=0)
  Mstar <- t(t(pstar)*Tdotk)
  uik <- pstar - censustab
  
  ##Generate sample design proportions
  MhatikMtildedes <- genmixeddirsnozerosspNEWbNRI(pstar, Tdotk, nks, ncs,deltasmat,alphasvec,taus,fakeck, K)
    
  genocslist.temp$phatik <- MhatikMtildedes[[1]]
  genocslist.temp$pstar <-pstar
  genocslist.temp$Pik <- censustab
  genocslist.temp$Sighateedir <- MhatikMtildedes[[2]]
  genocslist.temp$nks <- nks
  genocslist.temp$Tdotk <- Tdotk
  genocslist.temp
  
}



