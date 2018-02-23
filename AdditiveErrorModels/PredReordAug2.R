
##### Raking
predrakReord <- predfunrak( phatpredmod, datlistnew[[1]]$phatik, Tdotkreord)

##### Aug BHF
if(restrictthet == 1){
  predaugBHFReord  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunBHF, Tdotkreord)[[3]]    
}else{
  predaugBHFReord  <-  predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  bsiguhat[[1]][1,], datlistnew[[1]]$Sighateedir, psifunBHF, Tdotkreord)[[3]]    
  
}

##### Aug W-invserse
if(restrictthet == 1){
  predaugWinvReord  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunWinv, Tdotkreord)[[3]]    
}else{
  predaugWinvReord  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  bsiguhat[[1]][1,], datlistnew[[1]]$Sighateedir, psifunWinv, Tdotkreord)[[3]]    
}
 



#####  MSE BHF 
if(restrictthet == 1){
  MSEAugBHFReord <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotkreord, Xcatthet)[[1]])
}else{
  MSEAugBHFReord <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotkreord, Xcatthet)[[1]]) 
}

###### MSE Aug W-inverse
if(restrictthet == 1){
  MSEAugWinvReord <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotkreord, Xcatthet)[[1]])
}else{
  MSEAugWinvReord <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotkreord, Xcatthet)[[1]])
}

#####  MSE Raking
if(restrictthet == 1){
  MSERakReord <- diag(rakmse(datlistnew[[1]]$phatik, c(bsiguhat[[1]][1,1],1),  Tdotkreord, phatpredmod, kappahatcurrent,  Chatk[[1]], SighateemodfunLFS, nks, restrictthet, m, K, Xcatthet, datlistnew[[1]]$Sighateedir))
}else{
  MSERakReord <- diag(rakmse(datlistnew[[1]]$phatik, bsiguhat[[1]][1,],  Tdotkreord, phatpredmod, kappahatcurrent,  Chatk[[1]], SighateemodfunLFS, nks, restrictthet, m, K, Xcatthet, datlistnew[[1]]$Sighateedir))
}



