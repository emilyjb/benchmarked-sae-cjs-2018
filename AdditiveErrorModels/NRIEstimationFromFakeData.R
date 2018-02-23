
rm(list= ls(all = TRUE))
 

library("Matrix")
library("MASS")

load("FakeDataNRI.Rdata")


source("augfuns.R")
source("estfuns.R")
source("genfuns.R")


m <- 2; K <- 87


#####  Calculate variance model

pTinit.0 <- initbinom(datlistNRI[[1]]$phatik, datlistNRI[[1]]$Pik, datlistNRI[[1]]$nks, X, datlistNRI[[1]]$Tdotk,  startvals = rep(0, m), 1)
Guu <- gammauufun(pTinit.0[[1]])
vhatlist <- vector("list", K)
for(i in (1:K-1)*m+1){
	vhatlist[[i]] <- datlistNRI[[1]]$Sighateedir[i:(i+1),i:(i+1)]
}
VmatNRI <- matrix(unlist(vhatlist), nrow = 4, byrow = FALSE)
Ymatforvw <- t(apply(VmatNRI, 2, function(x){ diag(matrix(x, ncol = m)) }))[,-1]
Xmatforvw <- pTinit.0[[1]][1,]*pTinit.0[[1]][2,]/datlistNRI[[1]]$nks

chat.w <- lm(as.vector(as.matrix(Ymatforvw)-0.0001)~as.vector(Xmatforvw)-1)$coef
chat.w <- c(0.0001, chat.w)
c.vec <- rep(chat.w[2], length(datlistNRI[[1]]$nks))
datlistNRI[[1]]$Sighateedir <- SighateeFunNRI(chat.w, datlistNRI[[1]]$nks, pTinit.0[[1]])

##### Estimate parameters
pTinit <- initbinom(datlistNRI[[1]]$phatik, datlistNRI[[1]]$Pik, datlistNRI[[1]]$nks, X, datlistNRI[[1]]$Tdotk,  startvals = rep(0, m), 2)
betaEst <- estimatebetasigu11(datlistNRI, matrix(pTinit[[2]], nrow = 1), 0, 2,  chat.w, SighateeFunNRI, Inf, update.vardir = TRUE)
pFhat <- gfun(betaEst[[1]][1,], X, m, K)
datlistNRI[[1]]$Sighateedir <- SighateeFunNRI(chat.w, datlistNRI[[1]]$nks,pFhat)
varbetahat <- vhatbetahat(X, betaEst[[1]][1,], chat.w, SighateeFunNRI, datlistNRI[[1]]$nks, betaEst[[4]], 2, m, K, datlistNRI[[1]]$Sighateedir)
##### Standardized residuals
rstandard <- as.vector(datlistNRI[[1]]$phatik[2,] - pFhat[2,])/sqrt(betaEst[[4]] + matrix(diag(datlistNRI[[1]]$Sighateedir),m,K)[2,])
plot(matrix(X[,2],m,K)[2,], rstandard, xlab = expression(paste("Logit of CDL Proportion ",  x[ik])), ylab = "Standardized Residual")
abline(h = 0)

####Estimates of model parameters:
##### beta
betaEst[[1]]
sqrt(diag(varbetahat))
##### sigma2u
betaEst[[4]]
sqrt(betaEst[[5]])

##### Model-based (non-benchmarked predictors)
phatpredmod <- predfunmod(betaEst[[1]][1,], X, datlistNRI[[1]]$phatik, datlistNRI[[1]]$nks, betaEst[[4]], chat.w, SighateeFunNRI)
total.pred1 <- apply(t(t(phatpredmod)*datlistNRI[[1]]$Tdotk/sum(datlistNRI[[1]]$Tdotk)), 1, sum)
total.pred2 <- apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk/sum(datlistNRI[[1]]$Tdotk)), 1, sum)
total.pred1 - total.pred2
MSEhatpredmod <- vhatpredmod(X, betaEst[[1]][1,], chat.w, SighateeFunNRI, datlistNRI[[1]]$nks, betaEst[[4]],2, m, K, datlistNRI[[1]]$Sighateedir) 

##########    Evaluate weighted sum defining restriction
####################   Original model-based predictors (not benchmarked)
apply(t(t(phatpredmod)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk) 
####################   Direct estimators 
apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk) 



##### Aug-BHF predictor

predaugBHF <- predfunaug2(datlistNRI[[1]]$phatik ,  X, chat.w, SighateeFunNRI, betaEst[[4]],restrictthet=2, datlistNRI[[1]]$nks, as.vector(betaEst[[1]]) , datlistNRI[[1]]$Sighateedir, psifunBHF, datlistNRI[[1]]$Tdotk) 

############### Confirm that benchmarking restriction is satisfied
apply(t(t(predaugBHF[[3]])*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(predaugBHF[[3]])*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)

apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)

####### MSE of Aug-BHF predictor
MSEAugBHF <- MSEaugPSIfunModified23(datlistNRI[[1]]$phatik,  chat.w, SighateeFunNRI, betaEst[[4]], 2, datlistNRI[[1]]$nks,  psifunBHF, datlistNRI[[1]]$Sighateedir,   betaEst[[1]][1,] , datlistNRI[[1]]$Tdotk, X)
diag(MSEAugBHF[[1]])

##### Raking predictor

predrak <- predfunrak( phatpredmod, datlistNRI[[1]]$phatik, datlistNRI[[1]]$Tdotk)
MSERak <- rakmse(datlistNRI[[1]]$phatik, betaEst[[1]][1,], datlistNRI[[1]]$Tdotk,phatpredmod, betaEst[[4]],  chat.w, SighateeFunNRI, datlistNRI[[1]]$nks, 2, m, K, X, datlistNRI[[1]]$Sighateedir)
diag(MSERak)
############### Confirm that benchmarking restriction is satisfied
apply(t(t(predrak[[2]])*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(predrak[[2]])*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)

apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)

#### Linear additive BHF psi

predBHFadd <- linaddbenchBHFadderror(betaEst[[1]][1,], X, m,K,datlistNRI[[1]]$Tdotk,datlistNRI[[1]]$phatik,phatpredmod , chat.w, SighateeFunNRI, datlistNRI[[1]]$nks, betaEst[[4]], datlistNRI[[1]]$Sighateedir)
MSEBHFadd <- MSELinAdd(datlistNRI[[1]]$phatik,   betaEst[[4]] , 2,  datlistNRI[[1]]$Sighateedir, X, betaEst[[1]][1,], datlistNRI[[1]]$nks, SighateeFunNRI, chat.w, datlistNRI[[1]]$Tdotk)

############### Confirm that benchmarking restriction is satisfied
apply(t(t(predBHFadd)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(predBHFadd)*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)

apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 1, sum)/sum(datlistNRI[[1]]$Tdotk)
apply(t(t(datlistNRI[[1]]$phatik)*datlistNRI[[1]]$Tdotk), 2, sum)/sum(datlistNRI[[1]]$Tdotk)


library("xtable")
rbind( summary(MSEBHFadd/matrix(diag(MSEhatpredmod), m, K)[-1,]), 
       summary(matrix(diag(MSERak) , m, K)[-1,]/matrix(diag(MSEhatpredmod), m, K)[-1,]),
       summary(matrix(diag(MSEAugBHF[[1]]), m, K)[-1,]/matrix(diag(MSEhatpredmod), m, K)[-1,])
     )
      



xtable(rbind( summary(MSEBHFadd/matrix(diag(MSEhatpredmod), m, K)[-1,]), 
       summary(matrix(diag(MSERak) , m, K)[-1,]/matrix(diag(MSEhatpredmod), m, K)[-1,]),
       summary(matrix(diag(MSEAugBHF[[1]]), m, K)[-1,]/matrix(diag(MSEhatpredmod), m, K)[-1,])
))

summary(MSEBHFadd/diag(datlistNRI[[1]]$Sighateedir))
sum(MSEBHFadd >diag( datlistNRI[[1]]$Sighateedir))
sum(diag(MSERak) >diag( datlistNRI[[1]]$Sighateedir))
sum(diag(MSEAugBHF[[1]]) >diag( datlistNRI[[1]]$Sighateedir))
