
library("xtable")

############################   Output for Supplement 

maxcnt <- iter

###### Supplement table set 1:

diffmserak <- apply((predraks[1:maxcnt,] - pstars[1:maxcnt,] )^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugd <- apply((predaugBHFs[1:maxcnt,]  - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseauge <- apply((predaugWinvs[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)

df.diffmsecompare <- data.frame( Ve = apply(phats[1:maxcnt,] - pstars[1:maxcnt,],2, var), predmod =apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean), predrak = diffmserak,  predaugd = diffmseaugd, predauge = diffmseauge)
xtable(data.frame(rep(nks, each = 2),rep( Tdotk, each= 2)/sum(Tdotk), as.vector(Pik), df.diffmsecompare*10000) , digits = 3)

###### Supplement table set 2: 

diffmserakneword <- apply((predrakReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugdneword <- apply((predaugBHFReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugeneword <- apply((predaugWinvReords[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)

df.diffmsecompareset2 <- data.frame(Ve = apply((phats[1:maxcnt,] - pstars[1:maxcnt,]), 2, var),  predmod = apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2, 2, mean),  predrak = diffmserakneword,  predaugd = diffmseaugdneword, predauge = diffmseaugeneword)
data.frame(rep(nks, each = 2),rep( Tdotkreord, each= 2)/sum(Tdotk), as.vector(Pik), df.diffmsecompareset2*10000) 
xtable(data.frame(rep(nks, each = 2),rep( Tdotkreord, each= 2)/sum(Tdotk), as.vector(Pik), df.diffmsecompareset2*10000) )

###### Supplement table set 3: 

diffmserakroot <- apply((predrakroots[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugdroot <- apply((predaugBHFroots[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)
diffmseaugeroot <- apply((predaugWinvroots[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)-apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2,2,mean)

df.diffmsecompareset3 <- data.frame(Ve = apply((phats[1:maxcnt,] - pstars[1:maxcnt,]), 2, var),  predmod = apply((predunimods[1:maxcnt,] - pstars[1:maxcnt,])^2, 2, mean),  predrak = diffmserakroot,  predaugd = diffmseaugdroot, predauge = diffmseaugeroot)
data.frame(rep(nks, each = 2),rep( Tdotkroot, each= 2)/sum(Tdotkroot), as.vector(Pik), df.diffmsecompareset3*10000) 
xtable(data.frame(rep(nks, each = 2),rep( Tdotkroot, each= 2)/sum(Tdotkroot), as.vector(Pik), df.diffmsecompareset3*10000) )


