



linaddbenchBHFadderror <- function(params, X, m,K,Tdotk,phatik,pred.mat, vmodelparms, vmodelfun, nks, initkappahat , sighateedir){
  phatT <- gfun(params, X, m, K)
  ve <- diag(vmodelfun(vmodelparms, nks, phatT))
  gamma <- initkappahat*as.vector(phatT*(1-phatT))/(initkappahat*as.vector(phatT*(1-phatT)) + ve  )
  indouts <- ((1:K)-1)*m + 1
  phi.inv.vec <- (gamma*diag(sighateedir))[-indouts]
  W.matforbench <- kronecker(Tdotk, diag(m-1))  
  Llb <- solve(t(W.matforbench)%*%diag(phi.inv.vec)%*%W.matforbench)%*%t(W.matforbench)
  betahat.la <- solve(t(W.matforbench)%*%diag(phi.inv.vec)%*%W.matforbench)%*%t(W.matforbench)%*%(as.vector(phatik[-1,]) - as.vector(pred.mat[-1,]))
  pred.bench <- as.vector(pred.mat[-1,]) + diag(phi.inv.vec)%*%W.matforbench%*%betahat.la 
  rbind(1 - as.vector(pred.bench), as.vector(pred.bench))
}





MSELinAdd <- function(phatik,  initkappahat, restrictthet = 1,  Sighateedir, Xcatthet, params, nks, vmodelfun, vmodelparms, Tdotk){
  m <- dim(phatik)[1]; K <- dim(phatik)[2]
  #Xcatthet <- getXmats(allocslist[[twodind]][[2]])[[1]]
  X <- Xcatthet
  if(restrictthet==1){
    X <- Xcatthet[,-m]
    if(m==2){
      X <- matrix(X)
    }
  }else{
    X <- Xcatthet
  } 
  phatT <- gfun(params, X, m, K)
  Sighatee <- vmodelfun(vmodelparms, nks, phatT)
  ve <- diag(Sighatee)
  gamma <- initkappahat*as.vector(phatT*(1-phatT))/(initkappahat*as.vector(phatT*(1-phatT)) + ve  )
  indouts <- ((1:K)-1)*m + 1
  phi.inv.vec <- (gamma*diag(Sighateedir))[-indouts]
  W.matforbench <- kronecker(Tdotk, diag(m-1))  
  Llb <- solve(t(W.matforbench)%*%diag(phi.inv.vec)%*%W.matforbench)%*%t(W.matforbench)
  #  nikmat <- allocslist[[twodind]][[4]]
  #  nks <- apply(nikmat,2,sum)
  Xcat <- kronecker(rep(1,K), diag(m))[,-1]
  Xprov <- kronecker(diag(K),rep(1,m))
  #Xcen <- cbind(Xprov, Xlogitctab)
  guufun.list.b <- lapply(as.list(1:ncol(phatT)), function(x){(diag(phatT[,x]) - phatT[,x]%*%t(phatT[,x]))[-1,-1]})
  Guuinit <- as.matrix(bdiag(guufun.list.b))
  indouts <- ((1:K)-1)*m + 1
  #invGuuinit <- Guuinit
  #invGuuinit[-indouts,-indouts] <- solve(Guuinit[-indouts,-indouts])
  #invGuuinit[indouts,indouts] <- 0
  Xc <- X[-indouts,]
  Dgc <- Guuinit
  Vc <- initkappahat*Guuinit + Sighatee[-indouts,-indouts] 
  Lbeta <- solve(t(Xc)%*%Dgc%*%solve(Vc)%*%Dgc%*%Xc)%*%t(Xc)%*%Dgc 
  Hbeta <- Dgc%*%Xc%*%Lbeta
  Mw2 <- diag(phi.inv.vec)%*%W.matforbench%*%Llb
  matout <- Mw2%*%(Vc - Hbeta)%*%t(Mw2)
  
  #  mat1 <- diag(dim(Guuinit)[1]) - Dgc%*%Xc%*%Lbeta
  #  linmat <- Llb%*%diag(1 - gamma)%*%mat1
  #  matout <- linmat%*%Sighateedir[-indouts,-indouts]%*%t(linmat)
  term.lin <- diag(matout)
  term.mod <- vhatpredmod(Xcatthet, params, vmodelparms, vmodelfun, nks, initkappahat, restrictthet, m, K, Sighateedir)
  diag(term.mod)[-indouts] + term.lin
  
}

predfunmod <- function(params, Xcatthet, phatik, nks, kappahatcurrent, vmodelparms, vmodelfun){
  m <- dim(phatik)[1]; K <- dim(phatik)[2]
  phatT <- gfun(params, Xcatthet, m, K)
  ve <- diag(vmodelfun(vmodelparms, nks, phatT))
  gamma <- kappahatcurrent*as.vector(phatT*(1-phatT))/(kappahatcurrent*as.vector(phatT*(1-phatT)) + ve  )
  matrix(gamma*as.vector(phatik ) + (1-gamma)*as.vector(phatT), m, K, byrow = FALSE)
  
}


predfunrak <- function(phatpredmod, phatik, Wk){
  spree2(t(t(phatpredmod)*Wk), t(t(phatik)*Wk), 10)
}


psifunBHF <- function(kappahat, params, Xcatthet, Sighateedir, Tdotk, vmodelparms, vmodelfun, nks, m, K ){
  phatT <- gfun(params, Xcatthet, m, K)
  GammahatuuTpoolmods <- gammauufun(phatT)
  Dvv <- kappahat*diag(GammahatuuTpoolmods) + diag(Sighateedir)
  S2e <- vmodelfun(vmodelparms, nks, phatT)
  ve <- diag(S2e)
  gamma <- kappahat*as.vector(phatT*(1-phatT))/(kappahat*as.vector(phatT*(1-phatT)) + ve  )
  BHFphi.inv <- diag((1-gamma)^(-2))*kappahat*diag(GammahatuuTpoolmods)*diag(Sighateedir)/Dvv
  BHFphi.inv
  }

psifunWinv <- function(kappahat, params, Xcatthet, Sighateedir , Tdotk,vmodelparms, vmodelfun, nks, m, K){
  Constphi.inv <- (1/(rep(Tdotk,each = m)))
  phatT <- gfun(params, Xcatthet, m, K)
  GammahatuuTpoolmods <- gammauufun(phatT)
  S2e <- vmodelfun(vmodelparms, nks, phatT)
  ve <- diag(S2e)
  gamma <- kappahat*as.vector(phatT*(1-phatT))/(kappahat*as.vector(phatT*(1-phatT)) + ve  )
  Wphi.inv <- diag((1-gamma)^(-2)*Constphi.inv)
  Wphi.inv
}

gfunsmallmodX1X2b <- function(params1, params2, X1, X2, m, K){
  X <- cbind(X1, X2)
  params <- c(params1, params2)
  gfun(params, X, m, K)
}


rakmse <- function(phatik, params, Tdotk, predunimod,kappahatinits,  vmodelparms, vmodelfun,nks, restrictthet, m, K, Xcatthet, Sighateedir){
  ltab <- t(t(phatik)*Tdotk)
  ytotals <- as.vector(ltab)
  Mhatpredmod <- t(t(predunimod)*Tdotk)
  phirak <- thetatotals <- as.vector(Mhatpredmod)
  Xtilde1 <- kronecker(rep(1, K), diag(m))
  Xtilde2 <- kronecker(diag(K), rep(1, m))
  Xtilde <- cbind( Xtilde1, Xtilde2[,-1]   )#*Mhatrakmods[,iter]
  betahatlarak <- solve(t(Xtilde)%*%diag(phirak)%*%Xtilde)%*%t(Xtilde)%*%(ytotals - thetatotals)
  totalsrak1 <- diag(phirak)%*%Xtilde%*%betahatlarak%*%t(betahatlarak)%*%t(diag(phirak)%*%Xtilde)
  #msehat1mod <- MSEhatmod(Xcatthet, params, vmodelparms, nks, vmodelparms, vmodelfun, restrictthet)
  msehat1mod <- vhatpredmod(Xcatthet, params, vmodelparms, vmodelfun, nks, kappahatinits, restrictthet, m, K, Sighateedir)
  mserak1 <- msehat1mod + diag((Lfun2(ltab))%*%totalsrak1%*%(t(Lfun2(ltab))))
  mserak1
}

vhatpredmod <- function(Xcatthet, params, vmodelparms, vmodelfun, nks, kappahatinits, restrictthet, m, K, Sighateedir){
  pThatik <- gfun(params, Xcatthet, m, K)
  GammahatuuT <- gammauufun(pThatik)
  Lbeta <- linlamfun11(Xcatthet, pThatik, GammahatuuT, vmodelparms, vmodelfun, nks, kappahatinits, restrictthet)
  indouts <-seq(1, m*(K-1) + 1, by = m) 
  Ve <- vmodelfun(vmodelparms, nks, pThatik)
  Gammamat <- kappahatinits*GammahatuuT[-indouts,-indouts]%*%solve( kappahatinits*GammahatuuT[-indouts,-indouts] + Ve[-indouts,-indouts])
  if(restrictthet == 1){
    Xcatthet <- Xcatthet[,-m]
    if(m == 2){
      Xcatthet <- matrix(Xcatthet, ncol = 1)
    }
  }
  Upart <- Gammamat - diag(m*K)[-indouts,-indouts] + (diag(m*K)[-indouts,-indouts] - Gammamat)%*%GammahatuuT[-indouts,-indouts]%*%Xcatthet[-indouts,]%*%Lbeta
  Epart <- Upart + diag(m*K)[-indouts,-indouts]
  Ag1 <- kronecker(diag(1,K), rbind(-1, diag(m)[-1,]))
  Ag1 <- Ag1[,-indouts]
  MSE11 <- kappahatinits*Upart%*%GammahatuuT[-indouts,-indouts]%*%t(Upart)  + Epart%*%Sighateedir[-indouts,-indouts]%*%t(Epart)
  MSEhat <- Ag1%*%MSE11%*%t(Ag1)
  MSEhat
  
}

MSEhatmod <- function(Xcatthet, params, Ck, nks, vmodelparms, vmodelfun, restrictthet){
  phatT <- gfun(params, Xcatthet, m, K)
  ve <- diag(vmodelfun(vmodelparms, nks, phatT))
  gammaikmod <- kappahatcurrent*as.vector(phatT*(1-phatT))/(kappahatcurrent*as.vector(phatT*(1-phatT)) + ve  )
  if(restrictthet == 1){
Xrt1 <- Xcatthet[,-m]
  if(m == 2){
    Xrt1 <- as.matrix(Xrt1)
  }
  }else{
  Xrt1 <- Xcatthet
  }
  
GammahatuuT <- gammauufun(phatT)
ckvec <- rep(vmodelparms,each=m)
nksvec <- rep(nks,each=m)
msehat11mod <- diag(ckvec/nksvec)%*%diag(as.vector(gammaikmod))%*%GammahatuuT
Vmod <- diag(ckvec/nksvec+kappahatcurrent)%*%GammahatuuT
Vsynmod <- GammahatuuT%*%Xrt1%*%solve(t(Xrt1)%*%GammahatuuT%*%ginv(Vmod)%*%GammahatuuT%*%Xrt1)%*%t(Xrt1)%*%GammahatuuT
msehat1mod <- diag(msehat11mod) + (1-as.vector(gammaikmod))^2*diag(Vsynmod)

}


g3g4fun <- function(VhatkappahatNM, vhatchatks, ckvec, ntildevec, kappahatcurrent, m,K, Xcatthet, params){
  phatT <- gfun(params, Xcatthet, m, K)
  GammahatuuT <- gammauufun(phatT)
  ghat3ikmodnum <- (VhatkappahatNM*ckvec^2/ntildevec^2 + kappahatcurrent^2/ntildevec^2*rep(vhatchatks,each=m))*diag(GammahatuuT)
  ghat3ikmodden <- (kappahatcurrent + ckvec/ntildevec)^3
  ghat3ikmod <- ghat3ikmodnum/ghat3ikmodden
  ghat3ikmod
}




########  New augmented model function
predfunaug2 <- function(phatik ,  X, vmodelparms,  vmodelfun, initkappahat,restrictthet=1, nks,  params, Sighateedir, Psifun, Tdotk){
  #twodind <- 1
  m <- dim(phatik)[1]
  initlams <- rep(0,m-1)
  Xcat <- kronecker(rep(1,K), diag(m))[,-1]
  Dbench <- diag(rep(Tdotk, each = m))
  phatT <- gfun(params, X, m, K)
  S2e <- vmodelfun(vmodelparms, nks, phatT)
  ve <- diag(S2e)
  gamma <- initkappahat*as.vector(phatT*(1-phatT))/(initkappahat*as.vector(phatT*(1-phatT)) + ve  )
  Hgamma <- diag(1-gamma)
  if(is.null(Psifun)){
    Psi <- diag(diag(Hgamma)^(-2))
  }else{
    Psi <- Psifun(initkappahat, params, X, Sighateedir, Tdotk,  vmodelparms, vmodelfun, nks, m, K)
  }
  #Psi <- Thatdotkvec#/rep(nks,each=m)
  Xaug <- Psi%*%Hgamma%*%Dbench%*%Xcat  
  Ztilde <- Hgamma%*%Dbench%*%Xcat
  indouts <- ((1:K)-1)*m + 1
  Ztilde1 <- Ztilde[-indouts,]
  Psi1 <- Psi[-indouts,][,-indouts] 
  #Xaug <- diag(rep(Ck/nks,each=m))%*%Xcat
  Xprov <- kronecker(diag(K),rep(1,m))
  #Xcen <- cbind(Xprov, Xlogitctab)
  Guuinit <- gammauufun(phatT)
  Guuinit1 <- Guuinit[-indouts,][,-indouts]
  invGuuinit1 <- solve(Guuinit1)
  ltab <- t(t(phatik)*Tdotk)
  TPaug <- glsiternewEE.PSIder3(phatik,X[-indouts,],initkappahat,as.matrix(Xaug[-indouts,]),Hgamma,nks,startvals=initlams,params,Psi, invGuuinit1, Ztilde1, Tdotk, restrictthet)
  phatT <- gfun(params, X, m, K)
  ve <- diag(vmodelfun(vmodelparms, nks, phatT))
  gamma <- initkappahat*as.vector(phatT*(1-phatT))/(initkappahat*as.vector(phatT*(1-phatT)) + ve  )
  gammamat <- matrix(gamma, dim(ltab)[1], dim(ltab)[2], byrow = FALSE)
  pred <- gammamat*ltab + (1-gammamat)*TPaug[[1]]
  list(TPaug, pred, prop.table(pred, 2), gammamat)

}

glsstep2EE.PSIder3 <- function( params, X, piktab,phatik,PSI){
  #Sigaamod <- diag(rep(Ck/nks,each=m) + kappa)
  #Dwinv <- diag(1/PSI)#%*%solve(Sigaamod)
  #X <- solve(PSI)%*%X
  m <- dim(phatik)[1]
  K <- dim(phatik)[2]
  indouts <- ((1:K)-1)*m + 1
  PSI1 <- PSI[-indouts,][,-indouts]
  paramsnew <- params + solve(t(X)%*%PSI1%*%X)%*%t(X)%*%(as.vector(phatik[-1,])- as.vector(piktab[-1,]))
  paramsnew
}


glsiternewEE.PSIder3 <- function(phatik, Xoriginal, kappa,XforEst1,Hgamma,nks,startvals,params1,PSI,invGuuinit1, Ztilde, Tdotk, restrictthet=1){
  m <- dim(phatik)[1];
  #print(paste(restrictthet))
  params <- startvals
  cnt <- 0
  indouts <- ((1:K)-1)*m + 1
  repeat{
    cnt <- cnt + 1
    Z1 <- invGuuinit1%*%XforEst1
    prednum1 <- exp(Xoriginal%*%params1 + Z1%*%params)
    prednum1mat <- matrix(as.vector(prednum1), m-1, K, byrow = FALSE)
    predden <- 1 + apply(prednum1mat, 2, sum)
    piktab1 <- matrix(prednum1/predden, m-1, K, byrow = FALSE)
    piktab <- rbind(1/predden, piktab1)
    indouts <- 
    paramsnew <- glsstep2EE.PSIder3( params, Ztilde,piktab, phatik, PSI)
    diff <- max(abs(paramsnew - params))
    params <- paramsnew
    print(paste(cnt))
    if(diff<10^-10 || cnt>100){break}
  }
  Tiktab <- matrix(as.vector(piktab)*rep(Tdotk,each=m),m,K)
  ltab <- t(t(phatik)*Tdotk)
  score <- t((kronecker(diag(K), rep(1,m))))%*%diag(rep(Tdotk, each = m))%*%Hgamma%*%(as.vector(phatik) - as.vector(piktab))
  check <- apply(matrix(diag(Hgamma),m,K)*(ltab - Tiktab),1,sum)
  list(Tiktab,piktab,c(params1, params),diff,score,check) 
}



MSEaugPSIfunModified23 <- function(phatik,  vmodelparms,vmodelfun, initkappahat, restrictthet, nks,  Psifun, Sighateedir,params, Tdotk, Xcatthet){
  m <- dim(phatik)[1]; K <- dim(phatik)[2]
  #Xcatthet <- getXmats(allocslist[[twodind]][[2]])[[1]]
  X <- Xcatthet
  if(restrictthet==1){
    X <- Xcatthet[,-m]
    if(m==2){
      X <- matrix(X)
    }
  }else{
    X <- Xcatthet
  }
  Xcat <- kronecker(rep(1,K), diag(m))[,-1]
  Thatdotkvec <- rep(Tdotk,each=m)/sum(Tdotk)
  phatT <- gfun(params, Xcatthet, m, K)
  Dbench <- diag(Thatdotkvec)
  #Xaug <- Dbench%*%Xcat   
  Se <- vmodelfun(vmodelparms, nks, phatT)
  ve <- diag(Se)
  gamma <- initkappahat*as.vector(phatT*(1-phatT))/(initkappahat*as.vector(phatT*(1-phatT)) + ve  )
  Hgamma <- diag(1-gamma)
  if(is.null(Psifun) ){
    Psi <- diag(diag(Hgamma)^(-2))
  }else{
    Psi <- Psifun(initkappahat, params, Xcatthet, Sighateedir, Tdotk, vmodelparms, vmodelfun, nks,  m, K)
  }
  #Psi <- Thatdotkvec#/rep(nks,each=m)
  Xaug <- Psi%*%Hgamma%*%Dbench%*%Xcat    
  #Xaug <- diag(rep(Ck/nks,each=m))%*%Xcat
  Xprov <- kronecker(diag(K),rep(1,m))
  #Xcen <- cbind(Xprov, Xlogitctab)
  Guuinit <- gammauufun(phatT)
  indouts <- ((1:K)-1)*m + 1
  Guuinit1 <- Guuinit[-indouts,][,-indouts]
  invGuuinit1 <- solve(Guuinit1)
  Xc <- X[-indouts,]
  Dgc <- Guuinit[-indouts,-indouts]
  Vc <- (initkappahat*Guuinit + Se)[-indouts,][,-indouts]
  #Vc <- ((sqrt(Hgamma))%*%Guuinit%*%sqrt(Hgamma))[-indouts,-indouts]
  Lbeta <- solve(t(Xc)%*%Dgc%*%solve(Vc)%*%Dgc%*%Xc)%*%t(Xc)%*%Dgc%*%solve(Vc)
  W <- solve(Psi)%*%Xaug
  M <- solve(t(W[-indouts,])%*%Psi[-indouts,][,-indouts]%*%W[-indouts,])%*%t(W[-indouts,])
  Ldelta <- M%*%(diag(m*K)[-indouts,][,-indouts] - Guuinit[-indouts,][,-indouts]%*%X[-indouts,]%*%Lbeta)      
  Lbetadelta <- rbind(Lbeta, Ldelta)
  Z <- invGuuinit1%*%Xaug[-indouts,]
  Ag1 <- kronecker(diag(1,K), rbind(-1, diag(m)[-1,]))
  Ag1 <- Ag1[,-indouts]
  Dg1 <- Ag1%*%Dgc
  LpTaug <- Dg1%*%cbind(Xc, Z)%*%Lbetadelta
  Lhat <- Lfun2(t(t(phatik)*Tdotk))
  #Sighateedir <- Lhat%*%allocslist[[twodind]][[3]]%*%t(Lhat)
  #Sighateedir <- allocslist[[twodind]][[10]]
  #VpTaug <- LpTaug%*%(Guuinit[-indouts,-indouts]*initkappahat + Sighateedir[-indouts,-indouts])%*%t(LpTaug)
  Upart <- -Hgamma%*%Ag1 + Hgamma%*%LpTaug
  Epart <- (diag(m*K) - Hgamma)%*%Ag1 + Hgamma%*%LpTaug
  Laug2 <- Hgamma%*%Dg1%*%Z%*%Ldelta
  MSEaug2 <- Laug2%*%(initkappahat*Guuinit[-indouts,-indouts]+Sighateedir[-indouts,-indouts])%*%t(Laug2)
  MSEaug <- initkappahat*Upart%*%Guuinit[-indouts,-indouts]%*%t(Upart) + Epart%*%Sighateedir[-indouts,-indouts]%*%t(Epart)
  list(MSEaug, MSEaug2, X)
}




