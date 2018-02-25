rm(list = ls(all = TRUE))
source("augfuns.R")
source("estfuns.R")
source("genfuns.R")

library("Matrix")
library("MASS")


#########  Establish the fixed parameters for the simulation. 
source("fixedpars.R")
#########  Empty matrices for output
source("storeoutput.R")

datlistgen <- vector("list", 5000)

####### To reproduce the results in the paper, load this set of generated data for sigma^2_u = 0.005
#################### load("Sig2u005datlistgen.Rdata")
####### To reproduce the results in the supplement for sigma^2_u = 0.06, load this set of generated data: 
 load("Sig2u06datlistgen.Rdata")

### restrictthet = 1 restricts the slope to be 1 (as is the model for SPREE)
### restrictthet != 1 leaves the slope unrestricted, as in paper. 
restrictthet <- 2


iter <- 0

repeat{
  
  
iter <- iter + 1

##### Generate data
##### Comment the line with "datlistgen[[iter]]" out if using a set of existing generated data:
### datlistgen[[iter]] <- lapply(1:4, gendatafun, Pik, Tdotk, nks, ncs, deltasmat, alphasvec, taus, fakeck,K)

datlistnew <- datlistgen[[iter]]

##### Estimate variance model
Chatk <- compCkinit(datlistnew, 1)

for(i in 1:4){
  datlistnew[[i]]$Sighateedir <-  putzerosinVhatsPostDef2(datlistnew[[i]]$Sighateedir, Chatk[[1]], c(Chatk[[2]][i,1], 1), Xcatthet, datlistnew[[1]]$nks, m, K)
}

##### Estimate beta and sigma2u
if(restrictthet == 1){
  bsiguhat <- estimatebetasigu11(datlistnew,  Chatk[[2]]  ,0, restrictthet, Chatk[[1]], SighateemodfunLFS, 2)
}else{
  bsiguhat <- estimatebetasigu11(datlistnew,  cbind(Chatk[[2]],1)  ,0, restrictthet, Chatk[[1]], SighateemodfunLFS)
}
##### kappahatcurrent is sigma^2_u
kappahatcurrent <- bsiguhat[[4]]
vhatsig2u <- bsiguhat[[5]]
vhatchatk <- Chatk[[3]]
if(restrictthet ==1){
  varhatbetahat <- vhatbetahat(Xcatthet, c(bsiguhat[[1]][1,], 1), Chatk[[1]], SighateemodfunLFS, nks, kappahatcurrent, 1, m, K, datlistnew[[1]]$Sighateedir)
}else{
  varhatbetahat <- vhatbetahat(Xcatthet, bsiguhat[[1]][1,] , Chatk[[1]], SighateemodfunLFS, nks, kappahatcurrent, 1, m, K, datlistnew[[1]]$Sighateedir)
}
  
##### Estimated MMSE predictor -- not benchmakred
if(restrictthet == 1){
   phatpredmod <- predfunmod(c(bsiguhat[[1]][1,],1), Xcatthet, datlistnew[[1]]$phatik, nks, bsiguhat[[4]], Chatk[[1]],  SighateemodfunLFS)
}else{
   phatpredmod <- predfunmod(bsiguhat[[1]][1,], Xcatthet , datlistnew[[1]]$phatik, nks, bsiguhat[[4]], Chatk[[1]],  SighateemodfunLFS)  
}

##### Raking
predrak <- predfunrak( phatpredmod, datlistnew[[1]]$phatik, Tdotk)

##### Aug BHF
if(restrictthet == 1){
     predaugBHF <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunBHF, Tdotk)[[3]]    
}else{
     predaugBHF <-  predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=2, nks,  bsiguhat[[1]][1,], datlistnew[[1]]$Sighateedir, psifunBHF, Tdotk)[[3]]    
}

##### Aug W-invserse
if(restrictthet == 1){
  predaugWinv <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunWinv, Tdotk)[[3]]    
}else{
  predaugWinv <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks, bsiguhat[[1]][1,], datlistnew[[1]]$Sighateedir, psifunWinv, Tdotk)[[3]]    
}

##### MSE no benchmarking
if(restrictthet == 1){
  msehatpredmod  <- diag(vhatpredmod(Xcatthet, c(bsiguhat[[1]][1,], 1), Chatk[[1]], SighateemodfunLFS, nks, kappahatcurrent, 1, m, K, datlistnew[[1]]$Sighateedir))
}else{
  msehatpredmod  <- diag(vhatpredmod(Xcatthet, bsiguhat[[1]][1,], Chatk[[1]], SighateemodfunLFS, nks, kappahatcurrent, 2, m, K, datlistnew[[1]]$Sighateedir))
}


#####  MSE BHF 
if(restrictthet == 1){
  MSEAugBHF <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotk, Xcatthet)[[1]])
}else{
  MSEAugBHF <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotk, Xcatthet)[[1]]) 
}

###### MSE Aug W-inverse
if(restrictthet == 1){
  MSEAugWinv <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotk, Xcatthet)[[1]])
}else{
  MSEAugWinv <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotk, Xcatthet)[[1]])
}

#####  MSE Raking
if(restrictthet == 1){
  MSERak <- diag(rakmse(datlistnew[[1]]$phatik, c(bsiguhat[[1]][1,1],1),  Tdotk, phatpredmod, kappahatcurrent,  Chatk[[1]], SighateemodfunLFS, nks, restrictthet, m, K, Xcatthet, datlistnew[[1]]$Sighateedir))
}else{
  MSERak <- diag(rakmse(datlistnew[[1]]$phatik,  bsiguhat[[1]][1,] ,  Tdotk, phatpredmod, kappahatcurrent,  Chatk[[1]], SighateemodfunLFS, nks, restrictthet, m, K, Xcatthet, datlistnew[[1]]$Sighateedir))
}

#### Reorder weights
Tdotkreord <- sort(Tdotk, decreasing = TRUE)
#### Square root of weights
Tdotkroot <- sqrt(Tdotk)

source("PredReordAug2.R")
source("PredRootAug2.R")


##### G3 term in MSE
if(restrictthet == 1){
    ghat3ik <- g3g4fun(vhatsig2u, vhatchatk,  rep(Chatk[[1]],  each = m), rep(datlistnew[[1]]$nks, each = m), kappahatcurrent, m, K, Xcatthet, c(bsiguhat[[1]][1,], 1))
}else{
  ghat3ik <- g3g4fun(vhatsig2u, vhatchatk,  rep(Chatk[[1]],  each = m), rep(datlistnew[[1]]$nks, each = m), kappahatcurrent, m, K, Xcatthet, bsiguhat[[1]][1,])
}


phats <- rbind(phats, as.vector(datlistnew[[1]]$phatik))
vhatedirects <- rbind(vhatedirects, diag(datlistnew[[1]]$Sighateedir))
pstars <- rbind(pstars, as.vector(datlistnew[[1]]$pstar))
sig2uhats <- c(sig2uhats, kappahatcurrent)
betahats <- rbind(betahats,bsiguhat[[1]][1,] )
varhatbetahats <- rbind(varhatbetahats,varhatbetahat )
predunimods <- rbind(predunimods, as.vector(phatpredmod ))
msehatpredmods <- rbind( msehatpredmods,  msehatpredmod)

predraks <- rbind(predraks, as.vector(predrak[[2]]))
predaugBHFs <- rbind(predaugBHFs, as.vector(predaugBHF))
predaugWinvs <- rbind(predaugWinvs, as.vector(predaugWinv))
MSEAugWinvs <- rbind(MSEAugWinvs, MSEAugWinv)
MSEAugRaks <- rbind(MSEAugRaks, MSERak)
MSEAugBHFs <- rbind(MSEAugBHFs, MSEAugBHF)

predrakReords <- rbind(predrakReords, as.vector(predrakReord[[2]]))
predaugBHFReords <- rbind(predaugBHFReords, as.vector(predaugBHFReord))
predaugWinvReords <- rbind(predaugWinvReords, as.vector(predaugWinvReord))
MSEAugWinvReords <- rbind(MSEAugWinvReords, MSEAugWinvReord)
MSEAugRakReords <- rbind(MSEAugRakReords, MSERakReord)
MSEAugBHFReords <- rbind(MSEAugBHFReords, MSEAugBHFReord)
ghat3iks <- rbind(ghat3iks,  ghat3ik )

checkmargBHF <- rbind(checkmargBHF, apply(t(predaugBHF)*Tdotk, 2, sum))
checkmargWinv <- rbind(checkmargWinv, apply(t(predaugWinv)*Tdotk, 2, sum))
margdirect <- rbind(margdirect, apply(t(datlistnew[[1]]$phatik)*Tdotk, 2, sum))

predrakroots <- rbind(predrakroots, as.vector(predrakroot[[2]]))
predaugBHFroots <- rbind(predaugBHFroots, as.vector(predaugBHFroot))
predaugWinvroots <- rbind(predaugWinvroots, as.vector(predaugWinvroot))
MSEAugWinvroots <- rbind(MSEAugWinvroots, MSEAugWinvroot)
MSEAugRakroots <- rbind(MSEAugRakroots, MSERakroot)
MSEAugBHFroots <- rbind(MSEAugBHFroots, MSEAugBHFroot)

print(paste("iter", iter))
 

if(iter%%100 == 0){
print(paste("Outer Count", iter))
}

if(iter == 5000){break}


}

maxcnt <- iter - 1

###### Table 2 set 1:

diffmserak <- apply((predraks[1:maxcnt,] - pstars[1:maxcnt,] )^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugd <- apply((predaugBHFs[1:maxcnt,]  - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseauge <- apply((predaugWinvs[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)

df.diffmsecompare <- data.frame( Ve = apply(phats[1:maxcnt,] - pstars[1:maxcnt,],2, var), predmod =apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean), predrak = diffmserak,  predaugd = diffmseaugd, predauge = diffmseauge)
 
###### Table 2 set 2: 

diffmserakneword <- apply((predrakReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugdneword <- apply((predaugBHFReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugeneword <- apply((predaugWinvReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)

df.diffmsecompareset2 <- data.frame(Ve = apply((phats[1:maxcnt,] - pstars[1:maxcnt,]), 2, var),  predmod = apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2, 2, mean),  predrak = diffmserakneword,  predaugd = diffmseaugdneword, predauge = diffmseaugeneword)
 
newtab2 <- rbind(data.frame(rep(nks, each = 2),Tdotk = rep( Tdotkreord, each= 2)/sum(Tdotkreord), as.vector(Pik), df.diffmsecompare*10000)[c(1,5,15,19),],
data.frame(rep(nks, each = 2),Tdotk = rep( Tdotk, each= 2)/sum(Tdotk), as.vector(Pik),  df.diffmsecompareset2*10000)[c(1,5,15,19),])

library("xtable")

xtable(data.frame(Set = c("Set 1", "Set 1", "Set 1", "Set 1", "Set 2", "Set 2", "Set 2", "Set 2"), newtab2))


##### Empirical coverages for augmented model predictors:
crits <- rep(qt(0.975, df =  nks),each=m)

lcids <- predaugBHFs[1:maxcnt,] -  t(sqrt(t(MSEAugBHFs[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
ucids <- predaugBHFs[1:maxcnt,] +  t(sqrt(t(MSEAugBHFs[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
augdset1cov <- tapply(matrix(apply(lcids < pstars[1:maxcnt,] & pstars[1:maxcnt,] < ucids, 2, mean),m,K)[1,],nks,mean)

lciauges <- predaugWinvs[1:maxcnt,] - t(sqrt(t(MSEAugWinvs[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)  
uciauges <- predaugWinvs[1:maxcnt,] + t(sqrt(t(MSEAugWinvs[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits) 
augeset1cov <- tapply(matrix(apply(lciauges < pstars[1:maxcnt,] & pstars[1:maxcnt,] < uciauges, 2, mean),m,K)[1,],nks,mean)

lciraks <- predraks[1:maxcnt,] - t(sqrt(t(MSEAugRaks[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)  
uciraks <- predraks[1:maxcnt,] + t(sqrt(t(MSEAugRaks[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
rakset1cov <- tapply(matrix(apply(lciraks < pstars[1:maxcnt,] & pstars[1:maxcnt,] < uciraks, 2, mean),m,K)[1,],nks,mean)


lcidreords <- predaugBHFReords[1:maxcnt,] -  t(sqrt(t(MSEAugBHFReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
ucidreords <- predaugBHFReords[1:maxcnt,] +  t(sqrt(t(MSEAugBHFReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
augdset1covreord <- tapply(matrix(apply(lcidreords < pstars[1:maxcnt,] & pstars[1:maxcnt,] < ucidreords, 2, mean),m,K)[1,],nks,mean)

lciaugereords <- predaugWinvReords[1:maxcnt,] - t(sqrt(t(MSEAugWinvReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)  
uciaugereords <- predaugWinvReords[1:maxcnt,] + t(sqrt(t(MSEAugWinvReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits) 
augeset1covreord <- tapply(matrix(apply(lciaugereords < pstars[1:maxcnt,] & pstars[1:maxcnt,] < uciaugereords, 2, mean),m,K)[1,],nks,mean)

lcirakreords <- predrakReords[1:maxcnt,] - t(sqrt(t(MSEAugRakReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)  
ucirakreords <- predrakReords[1:maxcnt,] + t(sqrt(t(MSEAugRakReords[1:maxcnt,]+ghat3iks[1:maxcnt,]))*crits)
rakset1covreord <- tapply(matrix(apply(lcirakreords < pstars[1:maxcnt,] & pstars[1:maxcnt,] < ucirakreords, 2, mean),m,K)[1,],nks,mean)

#####  Output for table 3

xtable(rbind(augdset1cov, augdset1covreord,  augeset1cov, augeset1covreord, rakset1cov,   rakset1covreord))
