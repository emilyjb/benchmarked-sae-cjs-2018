
##### Raking
predrakroot <- predfunrak( phatpredmod, datlistnew[[1]]$phatik, Tdotkroot)

##### Aug BHF
if(restrictthet == 1){
  predaugBHFroot  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunBHF, Tdotkroot)[[3]]  
}else{
  predaugBHFroot  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=2, nks,   bsiguhat[[1]][1,] , datlistnew[[1]]$Sighateedir, psifunBHF, Tdotkroot)[[3]]
  
}

##### Aug W-invserse
if(restrictthet == 1){
  predaugWinvroot  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=1, nks,  c(bsiguhat[[1]][1,1],1), datlistnew[[1]]$Sighateedir, psifunWinv, Tdotkroot)[[3]]  
}else{
  predaugWinvroot  <- predfunaug2(datlistnew[[1]]$phatik ,  Xcatthet, Chatk[[1]], SighateemodfunLFS, kappahatcurrent,restrictthet=2, nks,   bsiguhat[[1]][1,] , datlistnew[[1]]$Sighateedir, psifunWinv, Tdotkroot)[[3]]
}


#####  MSE BHF 
if(restrictthet == 1){
  MSEAugBHFroot  <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotkroot, Xcatthet)[[1]])
}else{
  MSEAugBHFroot <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunBHF, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotkroot, Xcatthet)[[1]] )
}

###### MSE Aug W-inverse
if(restrictthet == 1){
  MSEAugWinvroot <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,  c(bsiguhat[[1]][1,1],1), Tdotkroot, Xcatthet)[[1]])
}else{
  MSEAugWinvroot <- diag(MSEaugPSIfunModified23(datlistnew[[1]]$phatik,  Chatk[[1]], SighateemodfunLFS, kappahatcurrent, restrictthet, nks,  psifunWinv, datlistnew[[1]]$Sighateedir,   bsiguhat[[1]][1,] , Tdotkroot, Xcatthet)[[1]])
}

#####  MSE Raking
if(restrictthet == 1){
  MSERakroot <- diag(rakmse(datlistnew[[1]]$phatik, c(bsiguhat[[1]][1,1],1),  Tdotkroot, phatpredmod, kappahatcurrent, Chatk[[1]], SighateemodfunLFS,nks,  restrictthet, 2, 10, Xcatthet, datlistnew[[1]]$Sighateedir))
}else{
  MSERakroot <- diag(rakmse(datlistnew[[1]]$phatik, bsiguhat[[1]][1,],  Tdotkroot, phatpredmod, kappahatcurrent, Chatk[[1]], SighateemodfunLFS, nks, restrictthet, 2, 10, Xcatthet, datlistnew[[1]]$Sighateedir))
}



