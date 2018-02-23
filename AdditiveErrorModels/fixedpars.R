
####  Type of weight (proportional to nk -- set 1 or inversely proportional to nk -- set 2)

Wtype <- "Proportional"

### Establish dimensions of two-way tables (notation different from paper
######  m = number of rows
######  K = number of areas = number of provinces
m <- 2; K <- 10

##### Construct covariates (in this case, representation of a table of census proporitons)
Pik <- matrix(0, m, 10)
Pik[1,] <- c(0.2, 0.75, 0.5, 0.2, 0.75, 0.75, 0.5, 0.2, 0.75, 0.5)
Pik[2,] <- 1 - Pik[1,]
###### Note: first row of pcikcen is gi of Table 2

###### Construct weights for benchmarking
###########  These have an interpretation as the column margins of the two-way table
Tdotk <- rep(c(2000,10000,35000,100000),each=3)[-c(4,7)]
Tdotkmat <- matrix(rep(Tdotk,each=m),m,K)
Wk <- Tdotk/sum(Tdotk)
########### Note: Wk is the weight (wi) in Table 2

Tik <-  Pik*Tdotkmat
TikPik <- list(Tik,Pik)

######  Sample sizes (corresponds to ni in table 2)
nks <- rep(c(16,30,60,204),each=3)[-c(4,7)]


#### Switch order of Tdotk for W inversely proportion to n (set 2)

if(Wtype == "Inverse"){
  
  Tdotk <- sort(Tdotk, decreasing = TRUE)
  Tik <- t(t(Pik)*Tdotk)
  TikPik <- list(Tik, Pik)
}

#####  Parameters for simulation configuration (Appendix B, Berg and Fuller, 2012)
deltasmat <- matrix(c(1, 1.25, 1.15, 1, 1.25, 1.15, 1, 1.25, 1.15, 1, 1.25, 1.15,
                      1, 1.55, 1.95, 1, 1.55, 1.95, 1, 1.55, 1.95, 1, 1.55, 1.95),2,12,byrow=TRUE)
deltasmat[deltasmat==1.95] <- 6
deltasmat[deltasmat==1.15] <- 6
deltasmat[,c(1,4,7,10)] <- 1

deltasmat <- deltasmat[,-c(4,7)]


#####  Parameters to generate mixture of Dirichlet to represent two-stage cluster (Appendix B, Berg and Fuller, 2012)
deltasmat <- matrix(c(1, 1.25, 1.15, 1, 1.25, 1.15, 1, 1.25, 1.15, 1, 1.25, 1.15,
                      1, 1.55, 1.95, 1, 1.55, 1.95, 1, 1.55, 1.95, 1, 1.55, 1.95),m-1,12,byrow=TRUE)
deltasmat[deltasmat==1.95] <- 6
deltasmat[deltasmat==1.15] <- 6
deltasmat[,c(1,4,7,10)] <- 1

deltasmat <- deltasmat[,-c(4,7)]

alphasvec <- c(1, 0.9, 0.8, 1, 0.90, 0.8, 1, 0.90, 0.80, 1, 0.90, 0.80)[-c(4,7)]
alphasvec[alphasvec==0.8] <- 0.9
#alphasvec[c(3,6,9,12)] <- 0.75
ncs <-  c(1, 2, 3, 1, 2, 3, 1, 2,3,1, 2,3)[-c(4,7)]
ncs[ncs==3] <- 2
taus <- c(1000,8,2,1000,8,2,1000,8,2,1000,8,2)[-c(4,7)]
taus[taus==2] <- 20
taus[taus==8] <- 4
taus[taus==1000] <- NA
ind <- 1:K
if(m>2){  #{deltasmat <-  matrix(deltasmat,1,K)}
  deffmat <- matrix(unlist(lapply(as.list(1:length(ind)), compd2funsp, Pik[,ind], deltasmat, alphasvec, ncs, taus)),3,length(ind))[-1,]
  deffmat[,is.na(taus)] <- 1
  fakeck <- apply(deffmat, 2, mean, na.rm=TRUE)
  fakeck[is.na(fakeck)] <- 1
  
}else{
  deltasmat <- matrix(deltasmat,1,10)
  fakeck <- rep(1,10)
}


#####  kappa is sigma^2_u in notation of manuscript
#####  Used kappa = 0.005 and kappa = 0.06
kappa <- 0.06



######  Set up matrics of row and column indicator variables
Fprov <- kronecker(diag(K), matrix(1,1,m))
Fcat <- kronecker(matrix(1,1,K),diag(m))
F <- rbind(Fcat, Fprov)
provfac <- factor(rep(1:K,each=m))
catfac <- factor(rep(1:m,times=K))


#####  Obtain interactions in census table (xi in Table 1)
cenint <- getcenint(Tik)

#####  Establish model matrices
Xcatcen <- cbind(t(Fcat),cenint)
Xmats <- getXmats(Tik)
Xcatthet <- Xmats[[1]][,-1]; Xcen <- Xmats[[2]][,-1]

#####  Covariance matrices of u and e
GammauuTtv <- gammauufun(TikPik[[2]])
nkvec <- rep(nks, each=m)
Dnkvec <- diag(nkvec)
Sigmaeetv <- solve(sqrt(Dnkvec))%*%GammauuTtv%*%solve(sqrt(Dnkvec))*(1-kappa)
Dprovtv <- diag(rep(Tdotk,each=m))
Sigaatv1 <- diag(rep(Tdotk^2,each=m)/nkvec*as.vector(Pik))
Sigaatv2 <- Dprovtv%*%solve(Dnkvec)%*%GammauuTtv%*%solve(Dnkvec)%*%Dprovtv
Sigaatv <- Sigaatv1 + Sigaatv2
Sigmauutv <- kappa*GammauuTtv
VMhatTdiffSRS <- (Dprovtv%*%Sigmauutv%*%Dprovtv + Sigaatv)
Tdotkmat2d <- matrix(rep(apply(TikPik[[1]],2,sum),each=m),m,K)
gammaikmodtv <- rep(kappa/(kappa + 1/nks),each=m)
