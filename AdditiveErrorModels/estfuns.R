

spree2 <- function(census, sampletab, steps){
	m <- dim(census)[1]; K <- dim(census)[2]
	margrow <- apply(sampletab, 1, sum)
	margcol <- apply(sampletab, 2, sum)
	start <- census
	cnt <- 0
	repeat{
		cnt <- cnt+1
		mtilde <- spreestep(start, margrow, margcol)
	start <- mtilde
	margnew1 <- apply(mtilde,1,sum)
	margnew2 <- apply(mtilde,2,sum)
	cval <- max(c(abs(margnew1-margrow), abs(margnew2-margcol)))
	if(cval<0.01|cnt==steps){break}
	}
	mtildemarg <- apply(mtilde, 2, sum)
	mtildemargmat <- matrix(rep(mtildemarg,each=m),m,K)
	ptilde <- mtilde/mtildemargmat
	list(mtilde,ptilde,cnt)
}

spreestep <- function(ests0, marg1, marg2){
	k <- dim(ests0)[[2]]
	m <- dim(ests0)[[1]]
	marg1mat <- matrix(rep(marg1, each=k), m, k, byrow=TRUE)
	marg2mat <- matrix(rep(marg2, each=m), m, k, byrow=FALSE)
	estmarg1 <- apply(ests0, 1, sum)
	estmarg1 <- matrix(rep(estmarg1, each=k), m, k, byrow=TRUE)
	update1 <- ests0/estmarg1*marg1mat
	estmarg2 <- apply(update1, 2, sum)
	estmarg2 <- matrix(rep(estmarg2, each=m), m, k, byrow=FALSE)
	update2 <- update1/estmarg2*marg2mat
	update2
}

gammauufun <- function(pik){
    m <- dim(pik)[1]; K <- dim(pik)[2]
    Pikmat <- matrix(rep(as.vector(pik), times = K), m * K, K)
    Gmat <- kronecker(diag(K), rep(1, m))
    PikG <- Gmat * Pikmat
    PPprime <- PikG %*% t(PikG)
    Gammauu <- diag(as.vector(pik)) - PPprime
    Gammauu
}

getinitvalsRT1nks <- function(table, sampletable, nks, Thatdotkmat){
   counttab <- t(t(sampletable)*nks)
   ThatPhat1 <- spree2(table, counttab, 50)
   alphai.init <- (log(ThatPhat1[[2]])[,1] - log(ThatPhat1[[2]])[1,1])[-1]
   list(ThatPhat1[[2]]*Thatdotkmat, ThatPhat1[[2]], alphai.init)
}


getcenint <- function(cTab){
  	as.vector(t( t(log(cTab)) - log(cTab)[1,]) - log(cTab)[,1] + log(cTab)[1,1])
}

getXmats <- function(ctab){
        m <- dim(ctab)[1]
        cenint <- getcenint(ctab)
        Xprovcatthet <- cbind(kronecker(rep(1,K),diag(m)), kronecker(diag(K),rep(1,m)),cenint)
        Xcatthet <- cbind(kronecker(rep(1,K),diag(m)), cenint)
        list(Xcatthet, Xprovcatthet)
}

gfun <- function(params, X, m, K){
  logitmu <- X%*%params
	explogitmu <- matrix(exp(logitmu),m,K)
	mu <- prop.table(explogitmu,2)
	mu
}


lfun <- function(i, k, spreetots){
  m <- dim(spreetots)[1]
  K <- dim(spreetots)[2]
  dKk <- basisfun(K, k)
  dmi <- basisfun(m, i)
  Mdotk <- apply(spreetots, 2, sum)[k]
  Tik <- spreetots[i,k]
  vm <- (Mdotk*dmi - Tik*rep(1,m))/Mdotk^2
  lik <- kronecker(matrix(dKk,1,K), matrix(vm,1,m))
  #list(dKk, vm, lik)
  lik
  
}


Lk2fun <- function(Thatk){
  Tdotk <- sum(Thatk)
  PhatTk <- Thatk/Tdotk
  m <- length(PhatTk)
  mat1 <- diag(rep(1/Tdotk,times=m))
  mat2 <- 1/Tdotk*matrix(rep(PhatTk,times=m),m,m)
  Lk <- mat1 - mat2
  Lk
}

Lfun2 <- function(That){
  Thatlist <- as.list(as.data.frame(That))
  Lkbigmat <- as.matrix(bdiag(lapply(Thatlist, Lk2fun)))
  Lkbigmat
}

getmodifiedcovmatPostDef1k <- function(k,Sigdir, Sigmod, nks){
  ind1 <- (k-1)*m + 1; ind2 <- k*m
  Sigdirk <- Sigdir[ind1:ind2,ind1:ind2]
  Sigmodk <- Sigmod[ind1:ind2,ind1:ind2]
  checkz <- diag(Sigdirk)==0
  if(sum(checkz)>0){
    indcheckz <- (1:m)[checkz]
    Vbarmodzz <- mean(diag(Sigmodk)[indcheckz])
    S <- 1/2/nks[k]^2
    W0 <- S/Vbarmodzz
    W <- min(W0, 1)
    W*Sigmodk + (1-W)*Sigdirk
  }else{
    Sigdirk
  }
}


putzerosinVhatsPostDef2 <- function(Sighateedir, Ck, params, X , nks, m, K){
  Phatnew <- gfun(params, X, m, K)
  GammahatuuTj <- gammauufun(Phatnew)
  cnvec <- rep(Ck,each=m)/rep(nks,each=m)
  Sighateemod <- diag(sqrt(cnvec))%*%GammahatuuTj%*%diag(sqrt(cnvec))
  modifiedSigee <- getmodifiedcovmatPostDef1(Sighateedir, Sighateemod, nkvec)
  modifiedSigee
}

getmodifiedcovmatPostDef1 <- function(V, Sighateemod, nks){
  listofVs <- lapply(as.list(1:K), getmodifiedcovmatPostDef1k, V, Sighateemod, nks)
  Vmod <- as.matrix(bdiag(listofVs))
  Vmod
}


updatebinom1step <- function(params, X, gammauustep, nks, phatik, piktab){
	m <- dim(phatik)[1]
      Dwn <- diag(rep(nks ,each=m))
      paramsnew <- params + solve(t(X)%*%gammauustep%*%Dwn%*%X)%*%t(X)%*%Dwn%*%(as.vector(phatik)- as.vector(piktab))
      paramsnew
}
	

initbinom <- function(phatik, Pik, nks, X, Thatdotk,  startvals = rep(0, m), restrictthet){
	if(sum(startvals == 0) == length(startvals)){
		startvals1 <- getinitvalsRT1nks(Pik, phatik, nks, rbind(Thatdotk, Thatdotk))
		if(restrictthet == 1){
			startvals <- startvals1[[3]]
		}else{
			startvals <- c(startvals1[[3]], 1)
		}
	}
	X1 <- X
	if(restrictthet == 1){
		if(m > 2){
			X <- X1[,-m]
		}else{
			X <- matrix(X1[,-m], m*K, 1)
		}
	}
	params <- startvals
	iter <- 0
	repeat{
		iter <- iter + 1
		if(restrictthet == 1){
			piktab <- gfun(c(params,1), X1, m, K)
		}else{
			piktab <- gfun(params, X, m, K)
		}
		gammauu <- gammauufun(piktab)	
		parUpdate <- updatebinom1step(params, X, gammauu, nks, phatik, piktab)
		check <- max(abs(parUpdate - params))
		params <- parUpdate
            if(check < 10^-5 || iter>30){break}
        }

   	if(restrictthet!=1){
                piktab <- gfun(params, X1, m, K)
        }
        if(restrictthet==1){
                piktab <- gfun(c(params,1),X1, m, K)
        }
	list(piktab, params)
}

compCkinit <- function(datlist, restrictthet){
	Cknums <- c()
	Ckdens <- c()
	betinitmat <- c()
	twodcnt <- 0
	repeat{
		twodcnt <- twodcnt + 1
		Sighateedir1t <- datlist[[twodcnt]]$Sighateedir
		nks <- datlist[[twodcnt]]$nks
		m <- dim(datlist[[twodcnt]]$phatik)[1]
		K <- dim(datlist[[twodcnt]]$phatik)[2]
		numck <- tapply(rep(nks,each=m)*diag(Sighateedir1t), rep(1:K, each=m), sum)
		X <- getXmats(datlist[[twodcnt]]$Pik)[[1]][,-1]
		PhatTmatPar <- initbinom(datlist[[twodcnt]]$phatik, datlist[[twodcnt]]$Pik, datlist[[twodcnt]]$nks, X, datlist[[twodcnt]]$Tdotk,  startvals = rep(0, m), restrictthet)
		betinitmat <- rbind(betinitmat, PhatTmatPar[[2]])
		Gut <- gammauufun(PhatTmatPar[[1]])
		denck <- tapply(diag(Gut), rep(1:K, each=m), sum)
		Cknums <- cbind(Cknums, numck)
		Ckdens <- cbind(Ckdens, denck)
		if(twodcnt==length(datlistnew)){break}
	}	
	Ck <- apply(Cknums,1,sum)/apply(Ckdens,1,sum)
	##### Estimated variance of estimator of c_k
	Ckmat <- matrix(rep(Ck,times=length(datlist)),K,length(datlist))
	Ckallocs <- Cknums/Ckdens
	vhatchatks <- apply(Ckdens*(Ckallocs - Ckmat)^2,1,sum)/apply(Ckdens,1,sum)/3
	if(sum(Ck==0)>0){
	vhatchatks[Ck == 0] <- Inf
	}
	wck <- 10000/(vhatchatks+ 10000)
	Ck <- Ck*wck + (1-wck)*1
	#vhatchatks <- vhatchatksall[,1]
	list(Ck, betinitmat, vhatchatks)
}



glsstep2EE <- function(params, X,gammauustep,kappa,vmodelparms, vmodelfun, nks,piktab,phatik){
  Vee <- vmodelfun(vmodelparms, nks, piktab)
  Vuu <- gammauustep*kappa
   Vinv <- ginv(Vee + Vuu)
  paramsnew <- params + solve(t(X)%*%gammauustep%*%Vinv%*%gammauustep%*%X)%*%t(X)%*%gammauustep%*%Vinv%*%(as.vector(phatik)- as.vector(piktab))
  paramsnew
}

glsiternewEEpostdef2 <- function(Pik, phatik, Tdotk, X, kappa, vmodelparms, vmodelfun, nks, startvals=rep(0,m),restrictthet=1){
        m <- dim(Pik)[1];
        if(sum(startvals == 0) == length(startvals)){
        		startvals1 <- getinitvalsRT1nks(Pik, phatik, nks, rbind(Tdotk, Tdotk))
		if(restrictthet == 1){
			startvals <- startvals1[[3]]
		}else{
			startvals <- c(startvals1[[3]], 1)
		}
	  }
        X1 <- X                 
        if(restrictthet==1){
                if(m>2){
                        X <- X1[,-m]
                }else{
                        X <- matrix(X1[,-m],m*K,1)
                }
        }
        params <- startvals
        cnt <- 0
        repeat{
                cnt <- cnt + 1
                if(restrictthet!=1){
                        piktab <- gfun(params, X, m, K)
                }
                if(restrictthet==1){
                        piktab <- gfun(c(params,1), X1, m, K)
                }
                gammauu <- gammauufun(piktab)   
                paramsnew <- glsstep2EE(params, X,gammauu,kappa,vmodelparms, vmodelfun, nks,piktab,phatik)
                diff <- max(abs(paramsnew - params))
                params <- paramsnew
                if(diff<10^-5 || cnt>30){break}
        }
        if(restrictthet!=1){
                piktab <- gfun(params, X1, m, K)
        }
        if(restrictthet==1){
                piktab <- gfun(c(params,1),X1, m, K)
        }
        Tiktab <- matrix(as.vector(piktab)*rep(Tdotk,each=m),m,K)
        list(Tiktab,piktab,params,diff)
}

SighateemodfunLFS <- function(Ck, nks, pThatik){
	  GammahatuuT <- gammauufun(pThatik)
        Ckvec <- rep(Ck, each=m)
        nksvec <- rep(nks,each=m)
        Sighateemod <- diag(sqrt(Ckvec/nksvec))%*%GammahatuuT%*%diag(sqrt(Ckvec/nksvec))
	  Sighateemod
}



SighateeFunNRI <- function(vmodelparms, nks,pThatik ){
		W.vec <- rep(nks, each = m)
		c.vec <- rep(vmodelparms[2], length(W.vec))
		GammahatuuT <- gammauufun(pThatik)
            V <- vmodelparms[1]*diag(length(W.vec)) + diag(sqrt(c.vec/W.vec))%*%GammahatuuT%*%diag(sqrt(c.vec/W.vec))
		V
}


linlamfun <- function(Xcatthet, pThatik, GammahatuuT,  Ck, vmodelfun, nks, kappahatinits, restrictthet){
	  Sighateemod <- vmodelfun(Ck, nks, pThatik)
        Sigmahatuumod <- kappahatinits*GammahatuuT
        VmodT <- Sigmahatuumod + Sighateemod
       if(m>2 & restrictthet==1){
                X <- Xcatthet[,-m]
        }
        if(m==2 & restrictthet==1){
                X <- matrix(Xcatthet[,-m],m*K,1)
        }
        if(restrictthet!=1){
                X <- Xcatthet
        }
        Linlam <- solve(t(X)%*%GammahatuuT%*%ginv(VmodT)%*%GammahatuuT%*%X)%*%t(X)%*%GammahatuuT%*%ginv(VmodT)
	  Linlam
}


getsigucomponentsfromGLSests <- function(phatik,  pThatik, GammahatuuT, Sighateedir0,nks,  vmodelparms, vmodelfun, kappahatinits, Xcatthet, restrictthet=1){
        m <- dim(phatik)[1]; K <- dim(phatik)[2]
        ##Compute GLSM estimators 
        #CENSUStab <- allocslist[[twodind]][[2]]
        #XcatXcen <- getXmats(CENSUStab)
        #Xcatthet <- XcatXcen[[1]]; Xcen <- XcatXcen[[2]][,-1]
        #GammahatuuT <- gammauufun(pThatik)
	  Sighateemod <- vmodelfun(vmodelparms, nks, pThatik)
        Sigmahatuumod <- kappahatinits*GammahatuuT
        VmodT <- Sigmahatuumod + Sighateemod
 	  Linlam <- linlamfun(Xcatthet, pThatik, GammahatuuT, vmodelparms,  vmodelfun, nks, kappahatinits, restrictthet)
        VmodT <- Sigmahatuumod + Sighateemod
        if(m>2 & restrictthet==1){
                X <- Xcatthet[,-m]
        }
        if(m==2 & restrictthet==1){
                X <- matrix(Xcatthet[,-m],m*K,1)
        }
        if(restrictthet!=1){
                X <- Xcatthet
        }
        LPT <- GammahatuuT%*%X%*%Linlam
        WPT <- diag(m*K) - LPT
        Siga <- WPT%*%Sighateemod%*%t(WPT)
        Sigb <- WPT%*%GammahatuuT%*%t(WPT)
        Siga1 <- WPT%*%Sighateedir0%*%t(WPT)
        sigb <- diag(Sigb)
        siga <- diag(Siga)
        ys <- (as.vector(phatik) - as.vector(pThatik))^2 - diag(Siga1)
        xs <- sigb
        Vhat2 <- solve(t(xs)%*%ginv(2*Siga^2)%*%xs)
        Vhat2diag <- solve(t(xs)%*%ginv(2*diag(diag(Siga))^2)%*%xs)
        xi <- 0.5*sqrt(Vhat2[1,1])
        xidiag <- 0.5*sqrt(Vhat2diag[1,1])
        Vhatnew <- 2*(Siga + kappahatinits*Sigb)^2
        Vhatnewdir <- 2*(Siga1 + kappahatinits*Sigb)^2
        psitildeden <- t(xs)%*%ginv(Vhatnew)%*%xs
        psitildenum <- (t(xs)%*%ginv(Vhatnew)%*%ys)
        numforvhatpooled <- t(xs)%*%ginv(Vhatnew)%*%Vhatnewdir%*%ginv(Vhatnew)%*%xs
        initks <- c(1/Vhat2[1,1], psitildenum,psitildeden, numforvhatpooled)
        initks
}



estimatebetasigu <- function(datlistnew, betinits, kappahatinits, restrictthet, vmodelparms, vmodelfun, maxiter = 2, update.vardir = FALSE){
	outeriter <- 0
	repeat{
		betcuralls <- c()
		sigustuff <- c()
		outeriter <- outeriter + 1
		for(i in 1:length(datlistnew)){
			Xcatthet <- getXmats(datlistnew[[i]]$Pik)[[1]][,-1]
			PTlamhatupdate <- glsiternewEEpostdef2(datlistnew[[i]]$Pik, datlistnew[[i]]$phatik, datlistnew[[i]]$Tdotk, Xcatthet,kappahatinits, vmodelparms, vmodelfun, datlistnew[[i]]$nks, startvals=betinits[i,],restrictthet=restrictthet)
			betcur <- PTlamhatupdate[[3]]
			GammahatuuT <- gammauufun(PTlamhatupdate[[2]])	
			if(update.vardir){
			 datlistnew[[i]]$Sighateedir <- vmodelfun(vmodelparms, datlistnew[[i]]$nks, PTlamhatupdate[[2]])
			}
			sigustuffiter <- getsigucomponentsfromGLSests(datlistnew[[i]]$phatik,  PTlamhatupdate[[2]], GammahatuuT, datlistnew[[i]]$Sighateedir, datlistnew[[i]]$nks, vmodelparms, vmodelfun, kappahatinits, Xcatthet, restrictthet = restrictthet)
			sigustuff <- rbind(sigustuff, sigustuffiter)
			betcuralls <- rbind(betcuralls, as.vector(betcur))
		}
		siguupdate <- sum(sigustuff[,2])/sum(sigustuff[,3])
		vhatsiguupdate <- sum(sigustuff[,4])/sum(sigustuff[,3])^2
		maxdiff <- max(c(abs(betcuralls - betinits), abs(siguupdate - kappahatinits)) )
		betinits <- betcuralls
		kappahatinits <- max(c(siguupdate, 1/sum(sigustuff[,1])))
		if(maxdiff < 10^-10 | outeriter > maxiter){break}
	}
	list(betinits, kappahatinits, 1/sum(sigustuff[,1]),  max(c( kappahatinits, 0.5*sqrt(1/sum(sigustuff[,1])))), vhatsiguupdate,   outeriter)
}




SighateeFunNRI <- function(vmodelparms, nks,pThatik ){
  W.vec <- rep(nks, each = m)
  c.vec <- rep(vmodelparms[2], length(W.vec))
  GammahatuuT <- gammauufun(pThatik)
  V <- vmodelparms[1]*kronecker(diag(K), matrix(c(1, -1, -1, 1), 2, 2) ) + diag(sqrt(c.vec/W.vec))%*%GammahatuuT%*%diag(sqrt(c.vec/W.vec))
  V
}



vhatbetahat <- function(Xcatthet, params, vmodelparms, vmodelfun, nks, kappahatinits, restrictthet, m, K, Sighateedir){
  pThatik <- gfun(params, Xcatthet, m, K)
  GammahatuuT <- gammauufun(pThatik)
  Lbeta <- linlamfun11(Xcatthet, pThatik, GammahatuuT, vmodelparms, vmodelfun, nks, kappahatinits, restrictthet)
  indouts <-seq(1, m*(K-1) + 1, by = m) 
  kappahatinits*Lbeta%*%GammahatuuT[-indouts,-indouts]%*%t(Lbeta)  + Lbeta%*%Sighateedir[-indouts,-indouts]%*%t(Lbeta)
  
}


linlamfun11 <- function(Xcatthet, pThatik, GammahatuuT,  Ck, vmodelfun, nks, kappahatinits, restrictthet){
  Sighateemod <- vmodelfun(Ck, nks, pThatik) 
  Sigmahatuumod <- kappahatinits*GammahatuuT
  rm.row <- seq(1, m*(K-1) + 1, by = m)
  VmodT <- (Sigmahatuumod + Sighateemod)[-rm.row,-rm.row]
  if(m>2 & restrictthet==1){
    X <- Xcatthet[,-m]
  }
  if(m==2 & restrictthet==1){
    X <- matrix(Xcatthet[,-m],m*K,1)
  }
  if(restrictthet!=1){
    X <- Xcatthet
  }
  X <- X[-rm.row,]
  GammahatuuT <- GammahatuuT[-rm.row,-rm.row]
  Linlam <- solve(t(X)%*%GammahatuuT%*%solve(VmodT)%*%GammahatuuT%*%X)%*%t(X)%*%GammahatuuT%*%solve(VmodT)
  Linlam
}

getsigucomponentsfromGLSests11 <- function(phatik,  pThatik, GammahatuuT, Sighateedir0,nks,  vmodelparms, vmodelfun, kappahatinits, Xcatthet, restrictthet=1){
  m <- dim(phatik)[1]; K <- dim(phatik)[2]
  ##Compute GLSM estimators 
  #CENSUStab <- allocslist[[twodind]][[2]]
  #XcatXcen <- getXmats(CENSUStab)
  #Xcatthet <- XcatXcen[[1]]; Xcen <- XcatXcen[[2]][,-1]
  #GammahatuuT <- gammauufun(pThatik)
  Sighateemod <- vmodelfun(vmodelparms, nks, pThatik)
  Sigmahatuumod <- kappahatinits*GammahatuuT
  VmodT <- Sigmahatuumod + Sighateemod
  rm.row <- seq(1, m*(K-1) + 1, by = m)
  Linlam <- linlamfun11(Xcatthet, pThatik, GammahatuuT, vmodelparms,  vmodelfun, nks, kappahatinits, restrictthet)
  VmodT <- (Sigmahatuumod + Sighateemod)[-rm.row,-rm.row]
  if(m>2 & restrictthet==1){
    X <- Xcatthet[,-m]
  }
  if(m==2 & restrictthet==1){
    X <- matrix(Xcatthet[,-m],m*K,1)
  }
  if(restrictthet!=1){
    X <- Xcatthet
  }
  X <- X[-rm.row,]
  LPT <- (GammahatuuT[-rm.row,-rm.row]%*%X%*%Linlam) 
  WPT <- diag(m*K)[-rm.row,-rm.row] - LPT
  Siga <- WPT%*%Sighateemod[-rm.row,-rm.row]%*%t(WPT)
  Sigb <- WPT%*%GammahatuuT[-rm.row,-rm.row]%*%t(WPT)
  Siga1 <- WPT%*%Sighateedir0[-rm.row,-rm.row]%*%t(WPT)
  sigb <- diag(Sigb)
  siga <- diag(Siga)
  ys <- (as.vector(phatik[-1,]) - as.vector(pThatik[-1,]))^2 - diag(Siga1)
  xs <- sigb
  Vhat2 <- solve(t(xs)%*%ginv(2*Siga^2)%*%xs)
  Vhat2diag <- solve(t(xs)%*%ginv(2*diag(diag(Siga))^2)%*%xs)
  xi <- 0.5*sqrt(Vhat2[1,1])
  xidiag <- 0.5*sqrt(Vhat2diag[1,1])
  Vhatnew <- 2*(Siga + kappahatinits*Sigb)^2
  Vhatnewdir <- 2*(Siga1 + kappahatinits*Sigb)^2
  psitildeden <- t(xs)%*%ginv(Vhatnew)%*%xs
  psitildenum <- (t(xs)%*%ginv(Vhatnew)%*%ys)
  numforvhatpooled <- t(xs)%*%ginv(Vhatnew)%*%Vhatnewdir%*%ginv(Vhatnew)%*%xs
  initks <- c(1/Vhat2[1,1], psitildenum,psitildeden, numforvhatpooled)
  initks
}



estimatebetasigu11 <- function(datlistnew, betinits, kappahatinits, restrictthet, vmodelparms, vmodelfun, maxiter = 2, update.vardir = FALSE){
  outeriter <- 0
  repeat{
    betcuralls <- c()
    sigustuff <- c()
    outeriter <- outeriter + 1
    for(i in 1:length(datlistnew)){
      Xcatthet <- getXmats(datlistnew[[i]]$Pik)[[1]][,-1]
      PTlamhatupdate <- glsiternewEEpostdef2(datlistnew[[i]]$Pik, datlistnew[[i]]$phatik, datlistnew[[i]]$Tdotk, Xcatthet,kappahatinits, vmodelparms, vmodelfun, datlistnew[[i]]$nks, startvals=betinits[i,],restrictthet=restrictthet)
      betcur <- PTlamhatupdate[[3]]
      GammahatuuT <- gammauufun(PTlamhatupdate[[2]])	
      if(update.vardir){
        datlistnew[[i]]$Sighateedir <- vmodelfun(vmodelparms, datlistnew[[i]]$nks, PTlamhatupdate[[2]])
      }
      sigustuffiter <- getsigucomponentsfromGLSests11(datlistnew[[i]]$phatik,  PTlamhatupdate[[2]], GammahatuuT, datlistnew[[i]]$Sighateedir, datlistnew[[i]]$nks, vmodelparms, vmodelfun, kappahatinits, Xcatthet, restrictthet = restrictthet)
      sigustuff <- rbind(sigustuff, sigustuffiter)
      betcuralls <- rbind(betcuralls, as.vector(betcur))
    }
    siguupdate <- sum(sigustuff[,2])/sum(sigustuff[,3])
    vhatsiguupdate <- sum(sigustuff[,4])/sum(sigustuff[,3])^2
    maxdiff <- max(c(abs(betcuralls - betinits), abs(siguupdate - kappahatinits)) )
    betinits <- betcuralls
    kappahatinits <- max(c(siguupdate, 1/sum(sigustuff[,1])))
    if(maxdiff < 10^-10 | outeriter > maxiter){break}
  }
  list(betinits, kappahatinits, 1/sum(sigustuff[,1]),  max(c( kappahatinits, 0.5*sqrt(1/sum(sigustuff[,1])))), vhatsiguupdate, outeriter)
}





 

