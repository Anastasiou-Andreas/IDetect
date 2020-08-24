install.packages("devtools")
library(devtools)

library(wbs)
library(breakfast)
library(changepoint)
library(changepoint.np)
library(cpm)
library(cumSeg)
library(FDRSeg)
library(l1tf)
library(not)
library(R.cache)
library(R.oo)
library(R.utils)
library(R.methodsS3)
library(stepR)
library(Rcpp)
library(segmented)
library(Segmentor3IsBack)
library(IDetect)

options(expressions = 500000)

inner_prod_cumsum <- function(x, y = cumsum(x)) {
  n <- length(x)
  sqrt( ( (n-1) : 1) / n / (1 : (n-1))) * y[1 : (n-1)] -
    sqrt( (1 : (n-1)) / n / ( (n-1):1)) * (y[n] - y[ 1 : (n-1)])
  
}

## Our simulations to compare our method with other methods in the literature

## This function is used to find the mean squared error of each detected result with the true signal.

mean.from.cpt <- function(x, cpt) {
  n <- length(x)
  len.cpt <- length(cpt)
  if (len.cpt) cpt <- sort(cpt)
  beg <- endd <- rep(0, len.cpt+1)
  beg[1] <- 1
  endd[len.cpt+1] <- n
  if (len.cpt) {
    beg[2:(len.cpt+1)] <- cpt+1
    endd[1:len.cpt] <- cpt
  }
  means <- rep(0, len.cpt+1)
  for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])
  rep(means, endd-beg+1)
}


rev.sim.small <- function(x, s3ib = 15, ave_run2 = 3000, Kmax_wbs = 50, qmax_NOT = 25) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric", mean="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  
  print("S3IB")
  
  z <- Segmentor(x, model=2, Kmax = s3ib)
  S3IB <- new("cpt.est")
  S3IB@nocpt <- SelectModel(z)-1
  S3IB@cpt <- as.numeric(z@breaks[S3IB@nocpt+1, 1:S3IB@nocpt])
  S3IB@mean <- mean.from.cpt(x, S3IB@cpt)
  if (S3IB@cpt == length(x)){S3IB@cpt=0
  S3IB@nocpt=0}
  S3IB@time <- system.time(Segmentor(x, model=2, Kmax = s3ib))[[3]]
  
  print("pelt")
   
  z <- cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT")
  pelt <- new("cpt.est")
  if (z@cpts[1:(length(z@cpts)-1)] == length(x)){pelt@cpt = 0
  pelt@nocpt = 0}else{
    pelt@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)])
    pelt@nocpt <- length(pelt@cpt)}
  pelt@mean <- mean.from.cpt(x, pelt@cpt)
  pelt@time <- system.time(cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT"))[[3]]
   
   
  print("pelt.np")
  
  z <- cpt.np(x/mad(diff(x)/sqrt(2)), method="PELT")
  pelt.np <- new("cpt.est")
  if (z@cpts[1:(length(z@cpts)-1)] == length(x)){pelt.np@cpt = 0
  pelt.np@nocpt = 0}else{
    pelt.np@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)])
    pelt.np@nocpt <- length(pelt.np@cpt)}
  pelt.np@mean <- mean.from.cpt(x, pelt.np@cpt)
  pelt.np@time <- system.time(cpt.np(x/mad(diff(x)/sqrt(2)), method="PELT"))[[3]]
  
  print("FDRSeg")
   
  z <- fdrseg(x)
  FDRSeg <- new("cpt.est")
  FDRSeg@nocpt <- length(z$left) - 1
  if(FDRSeg@nocpt == 0){FDRSeg@cpt = 0}
  else {FDRSeg@cpt= z$left[2:(FDRSeg@nocpt + 1)] - 1}
  FDRSeg@mean <- mean.from.cpt(x, FDRSeg@cpt)
  FDRSeg@time <- system.time(fdrseg(x))[[3]]
   
  print("CPM")
  z <- processStream(x,cpmType = "Student")
  CPM_points <- new("cpt.est")
  CPM_points@cpt <- as.numeric(z$changePoints)
  CPM_points@nocpt <- length(CPM_points@cpt)
  CPM_points@mean <- mean.from.cpt(x, CPM_points@cpt)
  CPM_points@time <- system.time(processStream(x,cpmType = "Student"))[[3]]
  
  print("CPM2")
  z <- processStream(x,cpmType = "Student",ARL0 = ave_run2)
  CPM_points2 <- new("cpt.est")
  CPM_points2@cpt <- as.numeric(z$changePoints)
  CPM_points2@nocpt <- length(CPM_points2@cpt)
  CPM_points2@mean <- mean.from.cpt(x, CPM_points2@cpt)
  CPM_points2@time <- system.time(processStream(x,cpmType = "Student",ARL0 = ave_run2))[[3]]
   
  print("cumSeg")
  
  z <- jumpoints(x)
  cumSeg <- new("cpt.est")
  if (z$n.psi==0){cumSeg@cpt==0
    cumSeg@nocpt=0}
  else{
    cumSeg@cpt <- as.numeric(z$psi)
    cumSeg@nocpt <- length(cumSeg@cpt)}
  cumSeg@mean <- mean.from.cpt(x, cumSeg@cpt)
  cumSeg@time <- system.time(jumpoints(x))[[3]]
  
  print("wbs")
  
  z <- wbs(x)
  cpt.z_1 = changepoints(z,th.const=1,Kmax=Kmax_wbs)
  wbsth10 <- new("cpt.est")
  if(is.na(cpt.z_1$cpt.th[[1]])){wbsth10@cpt = 0}
  else{wbsth10@cpt <- cpt.z_1$cpt.th[[1]]}
  wbsth10@nocpt <- cpt.z_1$no.cpt.th
  wbsth10@mean <- mean.from.cpt(x, wbsth10@cpt)
  wbsth10@time <- system.time(changepoints(wbs(x,5000),th.const=1,Kmax=Kmax_wbs))[[3]]
  
  wbssbic <- new("cpt.est")
  if(is.na(cpt.z_1$cpt.ic[[1]])) {wbssbic@cpt = 0}
  else {wbssbic@cpt= cpt.z_1$cpt.ic[[1]]}
  wbssbic@nocpt <- cpt.z_1$no.cpt.ic[[1]]
  wbssbic@mean <- mean.from.cpt(x, wbssbic@cpt) ## The computational time is the same as in wbs10
  
  print("not")
  z <- not(x,method="not",contrast="pcwsConstMean")
  cpt.ic = features(z,q.max=qmax_NOT)
  notic <- new("cpt.est")
  if(any(is.na(cpt.ic$cpt))) {notic@cpt = 0
  notic@nocpt <- 0}
  else {notic@cpt= cpt.ic$cpt
  notic@nocpt <- length(notic@cpt)}
  notic@mean <- mean.from.cpt(x, notic@cpt)
  notic@time <- system.time(features(not(x,method="not",contrast="pcwsConstMean"),q.max = qmax_NOT))[[3]]
  
  print("IsolateDetect")
  IsolateDetect <- new("cpt.est")
  z <-  ID_pcm(x)
  IsolateDetect@cpt <- as.numeric(z$cpt)
  IsolateDetect@nocpt <- z$no_cpt
  IsolateDetect@mean <- mean.from.cpt(x, IsolateDetect@cpt)
  IsolateDetect@time <- system.time(ID_pcm(x))[[3]]

 print("IsolateDetect_max_thr")
 IsolateDetect_max_thr <- new("cpt.est")
 z <-  win_pcm_th(x, thr_con = sqrt(3/2))
 IsolateDetect_max_thr@cpt <- as.numeric(z)
 IsolateDetect_max_thr@nocpt <- length(z)
 IsolateDetect_max_thr@mean <- mean.from.cpt(x, IsolateDetect_max_thr@cpt)
 IsolateDetect_max_thr@time <- system.time(win_pcm_th(x,thr_con = sqrt(3/2)))[[3]]
  
  
  print("IsolateDetect_sdll")
  IsolateDetect_sdll <- new("cpt.est")
  z <-  breakfast:::model.sdll(sol.idetect(x))
  IsolateDetect_sdll@cpt <- as.numeric(z$cpts)
  IsolateDetect_sdll@nocpt <- z$no.of.cpt
  IsolateDetect_sdll@mean <- mean.from.cpt(x, IsolateDetect_sdll@cpt)
  IsolateDetect_sdll@time <- system.time(model.sdll(sol.idetect(x)))[[3]]

  print("wbs2_sdll")
  wbs2_sdll <- new("cpt.est")
  z <-  model.sdll(sol.wbs2(x))
  wbs2_sdll@cpt <- as.numeric(z$cpts)
  wbs2_sdll@nocpt <- z$no.of.cpt
  wbs2_sdll@mean <- mean.from.cpt(x, wbs2_sdll@cpt)
  wbs2_sdll@time <- system.time(model.sdll(sol.wbs2(x)))[[3]]
  
  print("Bottom up")
  z1 = breakfast:::model.thresh(sol.tguh(x), th_const = 1)
  tguh <- new("cpt.est")
  if(length(z1$cpt)==0) {tguh@cpt = 0
  tguh@nocpt=0}
  else {tguh@cpt= z1$cpt
  tguh@nocpt <- z1$no.of.cpt}
  tguh@mean <- z1$est
  tguh@time <- system.time(breakfast:::model.thresh(sol.tguh(x), th_const = 1))[[3]]
  
  list(pelt=pelt, FDRSeg=FDRSeg, pelt.np=pelt.np, CPM_points=CPM_points,CPM_points2=CPM_points2,
       S3IB=S3IB, cumSeg=cumSeg, wbsth10=wbsth10, wbssbic=wbssbic,notic=notic, tguh=tguh,
       IsolateDetect = IsolateDetect, IsolateDetect_max_thr = IsolateDetect_max_thr, 
       IsolateDetect_sdll = IsolateDetect_sdll,wbs2_sdll = wbs2_sdll
    )
}

## This is the function which carried on the simulation study. It uses the previous function for 100 draws
## of different time series with the same prespecified signal. As a measure of the accuracy of the detected number,
## it gives as output the difference between the number of detected and the number of true change-points. 
## As a measure of the accuracy of the estimated locations, it gives as output the Hausdorff distance and the
## Mean Squared error between the true signal and the one created by the detections. This function also returns
## the computational time required by each method.

rev.sim.study.small <- function(signal, true.cpt=NULL,sigma, m = 100, seed = NULL, gen_s3ib=15,
                                gen_ave_run2 = 3000, gen_Kmax_wbs = 50, gen_qmax_NOT = 25) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  S3IB <- new("est.eval")
  pelt <- new("est.eval")
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  FDRSeg <- new("est.eval")
  cumSeg <- new("est.eval")
  wbsth10 <- new("est.eval")
  wbssbic <- new("est.eval")
  IsolateDetect <- new("est.eval")
  IsolateDetect_max_thr <- new("est.eval")
  IsolateDetect_sdll <- new("est.eval")
  wbs2_sdll <- new("est.eval")
  notic <- new("est.eval")
  tguh <- new("est.eval")
  
  
  no.of.cpt <- sum(abs(diff(signal)) > 0)
  n <- length(signal)
  ns <- max(diff(c(0,true.cpt,n)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    
    est <- rev.sim.small(x, s3ib = gen_s3ib, ave_run2 = gen_ave_run2, Kmax_wbs = gen_Kmax_wbs, qmax_NOT = gen_qmax_NOT)
    
    S3IB@dnc[i] <- est$S3IB@nocpt - no.of.cpt
    S3IB@mse[i] <- mean((est$S3IB@mean - signal)^2)
    S3IB@diff <- abs(matrix(est$S3IB@cpt,nrow=no.of.cpt,ncol=length(est$S3IB@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$S3IB@cpt),byr=F))
    S3IB@dh[i] <- max(apply(S3IB@diff,1,min),apply(S3IB@diff,2,min))/ns
    S3IB@time[i] <- est$S3IB@time
    
    pelt@dnc[i] <- est$pelt@nocpt - no.of.cpt
    pelt@mse[i] <- mean((est$pelt@mean - signal)^2)
    pelt@diff <- abs(matrix(est$pelt@cpt,nrow=no.of.cpt,ncol=length(est$pelt@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt@cpt),byr=F))
    pelt@dh[i] <- max(apply(pelt@diff,1,min),apply(pelt@diff,2,min))/ns
    pelt@time[i] <- est$pelt@time
     
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@mse[i] <- mean((est$CPM_points@mean - signal)^2)
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@mse[i] <- mean((est$CPM_points2@mean - signal)^2)
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    FDRSeg@dnc[i] <- est$FDRSeg@nocpt - no.of.cpt
    FDRSeg@mse[i] <- mean((est$FDRSeg@mean - signal)^2)
    FDRSeg@diff <- abs(matrix(est$FDRSeg@cpt,nrow=no.of.cpt,ncol=length(est$FDRSeg@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$FDRSeg@cpt),byr=F))
    FDRSeg@dh[i] <- max(apply(FDRSeg@diff,1,min),apply(FDRSeg@diff,2,min))/ns
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@mse[i] <- mean((est$pelt.np@mean - signal)^2)
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    cumSeg@dnc[i] <- est$cumSeg@nocpt - no.of.cpt
    cumSeg@mse[i] <- mean((est$cumSeg@mean - signal)^2)
    cumSeg@diff <- abs(matrix(est$cumSeg@cpt,nrow=no.of.cpt,ncol=length(est$cumSeg@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$cumSeg@cpt),byr=F))
    cumSeg@dh[i] <- max(apply(cumSeg@diff,1,min),apply(cumSeg@diff,2,min))/ns
    cumSeg@time[i] <- est$cumSeg@time
    
    wbsth10@dnc[i] <- est$wbsth10@nocpt - no.of.cpt
    wbsth10@mse[i] <- mean((est$wbsth10@mean - signal)^2)
    wbsth10@diff <- abs(matrix(est$wbsth10@cpt,nrow=no.of.cpt,ncol=length(est$wbsth10@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbsth10@cpt),byr=F))
    wbsth10@dh[i] <- max(apply(wbsth10@diff,1,min),apply(wbsth10@diff,2,min))/ns
    wbsth10@time[i] <- est$wbsth10@time
    
    wbssbic@dnc[i] <- est$wbssbic@nocpt - no.of.cpt
    wbssbic@mse[i] <- mean((est$wbssbic@mean - signal)^2)
    wbssbic@diff <- abs(matrix(est$wbssbic@cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=F))
    wbssbic@dh[i] <- max(apply(wbssbic@diff,1,min),apply(wbssbic@diff,2,min))/ns
    wbssbic@time[i] <- wbsth10@time[i]
    
    notic@dnc[i] <- est$notic@nocpt - no.of.cpt
    notic@mse[i] <- mean((est$notic@mean - signal)^2)
    notic@diff <- abs(matrix(est$notic@cpt,nrow=no.of.cpt,ncol=length(est$notic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$notic@cpt),byr=F))
    notic@dh[i] <- max(apply(notic@diff,1,min),apply(notic@diff,2,min))/ns
    notic@time[i] <- est$notic@time
    
    IsolateDetect@dnc[i] <- est$IsolateDetect@nocpt - no.of.cpt
    IsolateDetect@mse[i] <- mean((est$IsolateDetect@mean - signal)^2)
    IsolateDetect@diff <- abs(matrix(est$IsolateDetect@cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect@cpt),byr=F))
    IsolateDetect@dh[i] <- max(apply(IsolateDetect@diff,1,min),apply(IsolateDetect@diff,2,min))/ns
    IsolateDetect@time[i] <- est$IsolateDetect@time
    
    IsolateDetect_max_thr@dnc[i] <- est$IsolateDetect_max_thr@nocpt - no.of.cpt
    IsolateDetect_max_thr@mse[i] <- mean((est$IsolateDetect_max_thr@mean - signal)^2)
    IsolateDetect_max_thr@diff <- abs(matrix(est$IsolateDetect_max_thr@cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect_max_thr@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect_max_thr@cpt),byr=F))
    IsolateDetect_max_thr@dh[i] <- max(apply(IsolateDetect_max_thr@diff,1,min),apply(IsolateDetect_max_thr@diff,2,min))/ns
    IsolateDetect_max_thr@time[i] <- est$IsolateDetect_max_thr@time
    
    IsolateDetect_sdll@dnc[i] <- est$IsolateDetect_sdll@nocpt - no.of.cpt
    IsolateDetect_sdll@mse[i] <- mean((est$IsolateDetect_sdll@mean - signal)^2)
    IsolateDetect_sdll@diff <- abs(matrix(est$IsolateDetect_sdll@cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect_sdll@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$IsolateDetect_sdll@cpt),byr=F))
    IsolateDetect_sdll@dh[i] <- max(apply(IsolateDetect_sdll@diff,1,min),apply(IsolateDetect_sdll@diff,2,min))/ns
    IsolateDetect_sdll@time[i] <- est$IsolateDetect_sdll@time
    
    wbs2_sdll@dnc[i] <- est$wbs2_sdll@nocpt - no.of.cpt
    wbs2_sdll@mse[i] <- mean((est$wbs2_sdll@mean - signal)^2)
    wbs2_sdll@diff <- abs(matrix(est$wbs2_sdll@cpt,nrow=no.of.cpt,ncol=length(est$wbs2_sdll@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbs2_sdll@cpt),byr=F))
    wbs2_sdll@dh[i] <- max(apply(wbs2_sdll@diff,1,min),apply(wbs2_sdll@diff,2,min))/ns
    wbs2_sdll@time[i] <- est$wbs2_sdll@time
    
    tguh@dnc[i] <- est$tguh@nocpt - no.of.cpt
    tguh@mse[i] <- mean((est$tguh@mean - signal)^2)
    tguh@diff <- abs(matrix(est$tguh@cpt,nrow=no.of.cpt,ncol=length(est$tguh@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$tguh@cpt),byr=F))
    tguh@dh[i] <- max(apply(tguh@diff,1,min),apply(tguh@diff,2,min))/ns
    tguh@time[i] <- est$tguh@time
    
    gc()
  }
  
list(pelt=pelt,pelt.np=pelt.np,FDRSeg=FDRSeg, CPM_points=CPM_points,CPM_points2=CPM_points2,
     S3IB = S3IB,cumSeg=cumSeg, wbsth10=wbsth10, wbssbic=wbssbic,notic=notic,
     IsolateDetect_max_thr = IsolateDetect_max_thr, IsolateDetect = IsolateDetect,
     IsolateDetect_sdll = IsolateDetect_sdll, wbs2_sdll = wbs2_sdll,tguh=tguh)
}

seed.temp = 16
justnoise = rep(0,3000)
SIMR1.small = rev.sim.study.small(justnoise,sigma=1,true.cpt=c(0),seed=seed.temp)

blocks = c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40),rep(10.98,308),
           rep(-4.39,82),rep(3.29,430),rep(19.03,225),rep(7.68,41),rep(15.37,61),rep(0,389))

SIMR2.small = rev.sim.study.small(blocks,true.cpt = c(205,267,308,472,512,820,902,1332,1557,1598,1659),sigma=10,seed=seed.temp)

teeth10 = c(rep(0,11),rep(1,10),rep(0,10),rep(1,10),rep(0,10),rep(1,10),rep(0,10),rep(1,10),rep(0,10),rep(1,10),rep(0,10),rep(1,10),rep(0,10),rep(1,9))
SIMR3.small = rev.sim.study.small(teeth10,true.cpt=seq(11,131,10),0.4,seed=seed.temp)

stairs10 = c(rep(1,11),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,9))
SIMR4.small = rev.sim.study.small(stairs10,true.cpt=seq(11,141,10),sigma=0.3,seed=seed.temp)

middle <- c(rep(0, 1000), rep(1.5,20), rep(0,980))
SIMR5.small = rev.sim.study.small(middle,true.cpt=c(1000,1020),sigma=1,seed=seed.temp, gen_ave_run2 = 2000)

long_teeth = rep(c(rep(0,10),rep(3,10)),1000)
SIMR6.small = rev.sim.study.small(long_teeth,true.cpt=seq(10,19990,10),sigma=0.8,seed=seed.temp, m=100, gen_s3ib = 2000, gen_ave_run2 = 20000, gen_Kmax_wbs = 2000, gen_qmax_NOT = 2000)

long_stairs = numeric()
step = list()
for (i in 1:500){
  step[[i]]= rep(2*(i-1),20)
  long_stairs = c(long_stairs,step[[i]])
}
SIMR7.small = rev.sim.study.small(long_stairs,true.cpt =  seq(20,9980,20),sigma=1,seed=seed.temp, gen_s3ib = 500,
                                  gen_Kmax_wbs = 500, gen_qmax_NOT = 500, gen_ave_run2 = 10000)


## For piecewise-linear signals
library(earth)
library(stats)
library(IDetect)
library(freeknotsplines)
library(not)
devtools::install_github("hadley/l1tf")
library("l1tf")

get.signal <- function(model){
  
  if(model$cpt.type == "pcwsConstMean"){
    
    signal <- rep(0, model$n)
    
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    signal[segments[1,1]:segments[1,2]] <- model$start[1]
    
    for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
    
    
  }else if(model$cpt.type == "pcwsLinContMean"){
    
    signal <- rep(0, model$n)
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    
    slope <- model$start[2]
    signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
    
    for(j in 2:nrow(segments)) {
      
      slope <- slope +  model$jump.size[j-1]
      
      for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
    }
    
    
  } 
  return(signal)
}


sim.model <- function(model, sigma=1){
  get.signal(model) + sigma * rnorm(model$n)} 


## The following function is to get the estimated signal in the case of piecewise-linearity.
fit_lin_cont <- function(x, cpt) {
  lx = length(x)
  if (missing(cpt)) 
    cpt <- linear_find_changepoint_ic_fixedw3(x)
  if (!is.null(cpt)) 
    if (any(is.na(cpt))) 
      cpt <- cpt[!is.na(cpt)]
    cpt <- sort(unique(c(cpt, 0, lx)))
    fit <- rep(0, lx)
    cpt <- setdiff(cpt, c(0, lx))
    X <- bs(1:lx, knots = cpt, degree = 1, intercept = TRUE)
    fit <- lm.fit(X, x)$fitted.values
    return(fit)
}

### Functions to run NOT, trend-filtering and CPOP.

## For CPOP, one needs to first download the zip file from http://www.research.lancs.ac.uk/portal/en/datasets/cpop(56c07868-3fe9-4016-ad99-54439ec03b6c).html
## and save it to their C folder.
source("C:/CPOP_Rcode/Rcode/CPOP.R")
source("C:/CPOP_Rcode/Rcode/estimateCPOP.R")
dyn.load("C:/CPOP_Rcode/Rcode/prune2R.dll")
dyn.load("C:/CPOP_Rcode/Rcode/coeff.updateR.dll")

## run CPOP input data
CPOP.run = function(z, beta = NULL, sig = 1, useC = T) {
  if (length(beta) != 1) 
    beta = 2 * log(length(z))  ##BIC penalty
  ans <- CPOP(z, beta = beta, sigsquared = sig^2, useC = useC, useCprune = useC)
  CPS = ans[[2]]
  est = estimate.CPOP(z, CPS)
  
  return(list(cpt = CPS[-c(1, length(CPS))], fit = est$f))
}


## run trend-filtering input data
tf.run = function(z, sig = 1, lambdas = exp(seq(log(10), log(1000), length = 250))) {
  M = length(lambdas)
  n = length(z)
  SIC = rep(0, M)
  for (i in 1:M) {
    lambda <- lambdas[i]
    ans.tf <- l1tf(z, lambda)  ##l1tf not there
    SIC[i] = sum((z - ans.tf)^2)
    dans <- round(diff(ans.tf), digits = 5)
    CP <- c()
    for (ii in 1:(length(dans) - 1)) {
      if (dans[ii] != dans[ii + 1]) {
        CP <- c(CP, ii)
      }
    }
    SIC[i] = SIC[i]/sig^2 + log(n) * length(CP)
  }
  k = which.min(SIC)
  ans.tf <- l1tf(z, lambdas[k])  ##l1tf not there
  dans <- round(diff(ans.tf), digits = 5)
  CP <- c()
  for (ii in 1:(length(dans) - 1)) {
    if (dans[ii] != dans[ii + 1]) {
      CP <- c(CP, ii)
    }
  }
  
  return(list(cpt = CP, fit = ans.tf, lam = lambdas[k]))
  
}

## The SDLL related functions for ID
all.slopechanges.are.cpts <- function(x) {
  
  diff.x <- abs(diff(diff(x)))
  
  cpts <- which(diff.x > 0)
  no.of.cpt <- length(cpts)
  est <- x
  
  
  list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
  
  
}

idetect.th_linear <- function(x, sigma = stats::mad(diff(diff(x))) / sqrt(6), thr_const = 1.4,
                              thr_fin = sigma * thr_const * sqrt(2 * log(length(x))),
                              s = 1, e = length(x), points = 3, k_l = 1, k_r = 1) {
  l <- length(x)
  Res <- matrix(0, 1, 4)
  y <- c(0, cumsum(x))
  points <- as.integer(points)
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s < 2) {
    Res_fin <- matrix(0, 1, 4)
    cpt <- 0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- breakfast:::start_end_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- IDetect:::cumsum_lin(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        x_temp_l <- x[left_points[k_l]:e]
        ipcl <- IDetect:::cumsum_lin(x_temp_l)
        pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
        CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
        Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
        if (CUSUM_l[k_l] > thr_fin) {
          chp <- pos_l[k_l]
          indic <- 1
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- IDetect:::cumsum_lin(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          x_temp_l <- x[left_points[k_l]:e]
          ipcl <- IDetect:::cumsum_lin(x_temp_l)
          pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
          CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
          Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
            indic <- 1
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (indic == 1) {
        r <- idetect.th_linear(x, s = s, e = chp, points = points,
                               thr_fin = thr_fin, k_r = k_r, k_l = 1)
      } else {
        r <- idetect.th_linear(x, s = chp + 1, e = e, points = points,
                               thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
      }
      cpt <- c(chp, r[[1]])
      Res_fin <- rbind(Res, r[[2]])
    } else {
      cpt <- chp
      Res_fin <- Res
    }
  }
  cpt <- cpt[cpt != 0]
  Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
  return(list(changepoints = sort(cpt), full_information = Res_fin, y = y))
}

window.idetect.th_linear <- function(xd, sigma = stats::mad(diff(diff(xd))) / sqrt(6), thr_con = 1.4,
                                     c_win = 5000, w_points = 3) {
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= c_win) {
    u <- idetect.th_linear(x = xd, thr_const = thr_con, points = w_points)
    return(u)
  } else {
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- list()
    uaddition <- list()
    tsm[[1]] <- xd[1:c_win]
    ufin <- idetect.th_linear(tsm[[1]], thr_fin = t, points = w_points)
    uaddition[[1]] <- list()
    uaddition[[1]] <- idetect.th_linear(x = xd[(max(1, c_win - (10 * w_points) + 1)):min( (c_win + (10 * w_points)), lg)], thr_fin = t, points = 2)
    uaddition[[1]][[1]] <- uaddition[[1]][[1]] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,1] <- uaddition[[1]][[2]][,1] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,2] <- min(uaddition[[1]][[2]][,2] + c_win - (10 * w_points),min( (c_win + (10 * w_points)), lg))
    uaddition[[1]][[2]][,3] <- uaddition[[1]][[2]][,3] + c_win - (10 * w_points)
    ufin[[1]] <- c(ufin[[1]], uaddition[[1]][[1]])
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- list()
      u[[i]] <- idetect.th_linear(x = tsm[[i]], thr_fin = t, points = w_points)
      u[[i]][[1]] <- u[[i]][[1]] + (i - 1) * c_win
      u[[i]][[2]][,1] <- u[[i]][[2]][,1] + (i - 1) * c_win
      u[[i]][[2]][,2] <- min(u[[i]][[2]][,2] + (i - 1) * c_win, i * c_win)
      u[[i]][[2]][,3] <- u[[i]][[2]][,3] + (i - 1) * c_win
      uaddition[[i]] <- list()
      uaddition[[i]] <- idetect.th_linear(x = xd[(max(1, i * c_win - (10 * w_points) + 1)):(min(i * c_win + (10 * w_points), lg))], thr_fin = t, points = 2)
      uaddition[[i]][[1]] <- uaddition[[i]][[1]] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,1] <- uaddition[[i]][[2]][,1] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,2] <- min(uaddition[[i]][[2]][,2] + i * c_win - (10 * w_points), min(i * c_win + (10 * w_points), lg))
      uaddition[[i]][[2]][,3] <- uaddition[[i]][[2]][,3] + i * c_win - (10 * w_points)
      ufin[[1]] <- c(ufin[[1]],u[[i]][[1]], uaddition[[i]][[1]])
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- list()
    u[[K]] <- idetect.th_linear(tsm[[K]], thr_fin = t, points = w_points)
    u[[K]][[1]] <- u[[K]][[1]]  + (K - 1) * c_win
    u[[K]][[2]][,1] <- u[[K]][[2]][,1]  + (K - 1) * c_win
    u[[K]][[2]][,2] <- min(u[[K]][[2]][,2]  + (K - 1) * c_win,lg)
    u[[K]][[2]][,3] <- u[[K]][[2]][,3]  + (K - 1) * c_win
    ufin_cpt <- c(ufin[[1]], u[[K]][[1]])
    Res_fin <- matrix(0, 1, 4)
    Res_fin <- rbind(Res_fin, ufin[[2]], uaddition[[1]][[2]])
    if (K > 2){
      for (i in 2:(K-1)){Res_fin <- rbind(Res_fin,u[[i]][[2]], uaddition[[i]][[2]])}}
    Res_fin <- rbind(Res_fin, u[[K]][[2]])
    Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
    return(list(changepoints = sort(ufin_cpt), full_information = Res_fin,  y = c(0, cumsum(xd))))
  }
}

sol.idetect_linear <- function(x, thr_ic = 1.25, points = 3) {
  solutions.nested <- TRUE
  solution.set <- list()
  cands <- matrix(NA, 0, 4)
  lx <- length(x)
  if (lx < points) {solution.path <- integer()}
  else{
    points <- as.integer(points)
    step1 <- window.idetect.th_linear(x, thr_con = thr_ic, w_points = points)
    s1 <- as.matrix(step1$full_information)
    if (dim(s1)[2] == 1) {s1 <- t(s1)}
    ord <- order(s1[,4], decreasing = T)
    cands <- s1[ord, ,drop=FALSE]
    cpt_lower <- step1[[1]]
    lcpt_ic <- length(cpt_lower)
    seb_set <- c(unique(c(1, cpt_lower)), lx)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    CS <- matrix(cpt_lower,1,lcpt_ic)
    while (lseb_set >= 3) {
      Rs <- IDetect:::linear_contr_one(x, seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)],
                                       seb_set[2:(lseb_set - 1)])
      indic <- which.min(Rs)
      s1 <- setdiff(cpt_lower, seb_set[2:(lseb_set - 1)])
      d <- numeric(lcpt_ic)
      if(length(s1) == 0){d <- Rs}
      else{
        indic2 <- match(s1, cpt_lower)
        d[-indic2] <- Rs}
      CS <- rbind(CS,d)
      m_rs <- min(Rs)
      min_Rs <- seb_set[2:(lseb_set - 1)][indic]
      cands <- rbind(cands, c(seb_set[indic], seb_set[indic + 2], min_Rs, m_rs))
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    solution.path <- min_C[length(min_C):1]
    cusum_m <- apply(CS[-1,,drop = FALSE],2,max)
    indic3 <- match(cpt_lower, cands[,3])
    cands[indic3,4] <- cusum_m
    ord <- order(cands[,4], decreasing = T)
    cands <- cands[ord, ,drop=FALSE]#[-(length(solution.path)+1), ,drop = FALSE]
    cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    if(is.na(solution.path[1])){solution.path <- integer(0)}}
  ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, cands = cands, method = "idetect")
  
  class(ret) <- "cptpath"
  
  ret
}


model.sdll_linear <- function(cptpath.object, sigma = stats::mad(diff(diff(cptpath.object$x))) / sqrt(6), th.const = 1.25, th.const.min.mult = 0.3) {
  
  x <- cptpath.object$x
  
  n <- length(x)
  
  if (n <= 1) {
    
    est <- x
    
    no.of.cpt <- 0
    
    cpts <- integer(0)
    
  }
  
  else {
    
    if (sigma == 0) {
      
      s0 <- all.slopechanges.are.cpts(x)
      est <- s0$est
      no.of.cpt <- s0$no.of.cpt
      cpts <- s0$cpts
      
    } else {
      
        if (is.null(th.const)) {stop("th.const must be specified.")
        }
      }
      
      
      th.const.min <- th.const * th.const.min.mult
      
      th <- th.const * sqrt(2 * log(n)) * sigma
      
      th.min <- th.const.min * sqrt(2 * log(n)) * sigma
      
      
      if (cptpath.object$cands[1,4] < th) {
        
        no.of.cpt <- 0
        
        cpts <- integer(0)
        
      }
      
      else {
        
        indices <- which(cptpath.object$cands[,4] > th.min)
        
        if (length(indices) == 1) {
          
          cpts <- cptpath.object$cands[indices, 3]
          
          no.of.cpt <- 1
          
        }
        
        else {
          
          rc.sel <- cptpath.object$cands[indices,,drop=F]
          
          z <- cptpath.object$cands[indices,4]
          
          z.l <- length(z)
          
          dif <- -diff(log(z))
          
          dif.ord <- order(dif, decreasing=T)
          
          j <- 1
          
          while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1
          
          if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l
          
          cpts <- sort(cptpath.object$cands[1:no.of.cpt,3])			
          
        }
        
      } 
      
      est <- IDetect:::est_signal(x, cpts, type = "slope")
      
    }
  
  ret <- list(est=est, no.of.cpt=no.of.cpt, cpts=cpts, model="sdll", solution.path=cptpath.object$method)
  
  class(ret) <- "cptmodel"
  
  ret
  
}

## Simulation funtions

linear.rev.sim <- function(x, q_max_NOT = 25, FKS_knot = 10) {
  
  setClass("cpt.est", representation(cpt = "numeric", nocpt = "numeric", fit = "numeric", time = "numeric", dh = "numeric"), 
           prototype(cpt = numeric(0), nocpt = 0, fit = numeric(0), time = numeric(0), dh = numeric(0)))
  
  
  print("not")
  z <- not(x, method = "not", contrast = "pcwsLinContMean")
  cpt.ic = features(z, q.max = q_max_NOT)
  
  notic <- new("cpt.est")
  if (any(is.na(cpt.ic$cpt))) {
    notic@cpt = 0
    notic@nocpt = 0
  } else {
    notic@cpt = cpt.ic$cpt
    notic@nocpt <- length(notic@cpt)
  }
  notic@fit <- fit_lin_cont(x, notic@cpt)
  notic@time <- system.time(features(not(x, method = "not", contrast = "pcwsLinContMean"), q.max = q_max_NOT))[[3]]
  
  print("CPOP")
  z <- CPOP.run(x, sig = mad(diff(diff(x)))/sqrt(6), useC = F)

  cpop <- new("cpt.est")
  if (length(z$cpt) == 0) {
      cpop@cpt = 0
      cpop@nocpt = 0
  } else {
      cpop@cpt = z$cpt
      cpop@nocpt <- length(z$cpt)
  }
  cpop@fit <- z$fit
  cpop@time <- system.time(CPOP.run(x, sig = mad(diff(diff(x)))/sqrt(6), useC = F))[[3]]
  
  print("t1f")
  z <- tf.run(x/(mad(diff(diff(x)))/sqrt(6)))
  
  t1f <- new("cpt.est")
  if (length(z$cpt) == 0) {
      t1f@cpt = 0
      t1f@nocpt = 0
  } else {
      t1f@cpt = z$cpt
      t1f@nocpt <- length(z$cpt)
 }
  t1f@fit <- z$fit
  t1f@time <- system.time(tf.run(x/(mad(diff(diff(x)))/sqrt(6))))[[3]]
  
  print("MARS")
  z1 <- earth(1:length(x),y=x)
  z <- sort(unique(z1$cuts[z1$selected.terms,]))[-1]
  MARS <- new("cpt.est")
  if (length(z) == 0) {
    MARS@cpt = 0
    MARS@nocpt = 0
  }
  else {
    MARS@cpt = z
    MARS@nocpt = length(z)
  }
  MARS@fit <- fit_lin_cont(x, MARS@cpt)
  MARS@time <- system.time(sort(unique(earth(1:length(x),y=x)$cuts[earth(1:length(x),y=x)$selected.terms,]))[-1])[[3]]
  
  print("FKS")
  z1 <-  fit.search.numknots(1:length(x), y = x, degree = 1, minknot = 1, maxknot = FKS_knot,
                            alg = "LS", search = "genetic",
                            knotnumcrit = "adjGCV", k = 2, d = 3, seed = 5,
                            stream = 0)
  z <- round(z1@optknot)
  FKS <- new("cpt.est")
  if (length(z) == 0) {
    FKS@cpt = 0
    FKS@nocpt = 0
  } else {
    FKS@cpt = z
    FKS@nocpt <- length(z)
  }
  FKS@fit <- fit_lin_cont(x, FKS@cpt)
  FKS@time <- system.time(round(fit.search.numknots(1:length(x), y = x, degree = 1, minknot = 1, maxknot = FKS_knot,
                                                        alg = "LS", search = "genetic",
                                                        knotnumcrit = "adjGCV", k = 2, d = 3, seed = 5,
                                                        stream = 0)@optknot))[[3]]
  
   print("ID")
   IsolateDetect <- new("cpt.est")
   z = ID_cplm(x)
   IsolateDetect@cpt <- as.numeric(z$cpt)
   IsolateDetect@nocpt <- z$no_cpt
   IsolateDetect@fit <- fit_lin_cont(x, IsolateDetect@cpt)
   IsolateDetect@time <- system.time(ID_cplm(x))[[3]]
   
   print("ID_SDLL")
   IsolateDetect_sdll <- new("cpt.est")
   z = model.sdll_linear(sol.idetect_linear(x))
   IsolateDetect_sdll@cpt <- as.numeric(z$cpts)
   IsolateDetect_sdll@nocpt <- z$no.of.cpt
   IsolateDetect_sdll@fit <- IDetect:::est_signal(x, IsolateDetect_sdll@cpt, type = "slope")
   IsolateDetect_sdll@time <- system.time(model.sdll_linear(sol.idetect_linear(x)))[[3]]
  
  
  list(notic = notic, MARS = MARS, t1f = t1f, cpop = cpop, FKS = FKS, 
       IsolateDetect = IsolateDetect, IsolateDetect_sdll = IsolateDetect_sdll)
}


## This is the function which carried on the simulation study. It uses the previous function for 100 draws
## of different time series with the same prespecified signal. It gives as output the difference between
## the number of detected and the number of true change-points. It also gives as an output the Hausdorff distance
## and the Mean Squared error between the true signal and the one created by the detections.

linear_rev.sim.study <- function(model, sigma, m = 100, seed = NULL, gen_qmax = 25) {
  
  setClass("est.eval", representation(avg.signal = "numeric", fit = "list", cpt = "list", diff = "matrix", dh = "numeric", 
                                      cptall = "numeric", dnc = "numeric", mse = "numeric", time = "numeric"), prototype(dnc = numeric(m), mse = numeric(m), 
                                                                                                                         time = numeric(m), dh = numeric(m)))
  
  IsolateDetect <- new("est.eval")
  IsolateDetect_sdll <- new("est.eval")
  notic <- new("est.eval")
  cpop <- new("est.eval")
  t1f <- new("est.eval")
  MARS <- new("est.eval")
  FKS <- new("est.eval")
  
  signal = get.signal(model)
  
  no.of.cpt <- length(model$cpt)
  n <- length(signal)
  ns <- max(diff(c(0, model$cpt, n)))
  
  if (!is.null(seed)) 
    set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    
    est <- linear.rev.sim(x, q_max_NOT = gen_qmax, FKS_knot = no.of.cpt + 2)
    
    notic@dnc[i] <- est$notic@nocpt - no.of.cpt
    notic@mse[i] <- mean((signal - est$notic@fit)^2)
    notic@diff <- abs(matrix(est$notic@cpt, nrow = no.of.cpt, ncol = length(est$notic@cpt), byr = T) - matrix(model$cpt, 
                                                                                                              nrow = no.of.cpt, ncol = length(est$notic@cpt), byr = F))
    notic@dh[i] <- max(apply(notic@diff, 1, min), apply(notic@diff, 2, min))/ns
    notic@time[i] <- est$notic@time
    
    cpop@dnc[i] <- est$cpop@nocpt - no.of.cpt
    cpop@mse[i] <- mean((signal - est$cpop@fit)^2)
    cpop@diff <- abs(matrix(est$cpop@cpt, nrow = no.of.cpt, ncol = length(est$cpop@cpt), byr = T) - matrix(model$cpt, 
        nrow = no.of.cpt, ncol = length(est$cpop@cpt), byr = F))
    cpop@dh[i] <- max(apply(cpop@diff, 1, min), apply(cpop@diff, 2, min))/ns
    cpop@time[i] <- est$cpop@time
    
    t1f@dnc[i] <- est$t1f@nocpt - no.of.cpt
    t1f@mse[i] <- mean((signal - est$t1f@fit)^2)
    t1f@diff <- abs(matrix(est$t1f@cpt, nrow = no.of.cpt, ncol = length(est$t1f@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
        ncol = length(est$t1f@cpt), byr = F))
    t1f@dh[i] <- max(apply(t1f@diff, 1, min), apply(t1f@diff, 2, min))/ns
    t1f@time[i] <- est$t1f@time
    
    MARS@dnc[i] <- est$MARS@nocpt - no.of.cpt
    MARS@mse[i] <- mean((signal - est$MARS@fit)^2)
    MARS@diff <- abs(matrix(est$MARS@cpt, nrow = no.of.cpt, ncol = length(est$MARS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                           ncol = length(est$MARS@cpt), byr = F))
    MARS@dh[i] <- max(apply(MARS@diff, 1, min), apply(MARS@diff, 2, min))/ns
    MARS@time[i] <- est$MARS@time
    
    FKS@dnc[i] <- est$FKS@nocpt - no.of.cpt
    FKS@mse[i] <- mean((signal - est$FKS@fit)^2)
    FKS@diff <- abs(matrix(est$FKS@cpt, nrow = no.of.cpt, ncol = length(est$FKS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                        ncol = length(est$FKS@cpt), byr = F))
    FKS@dh[i] <- max(apply(FKS@diff, 1, min), apply(FKS@diff, 2, min))/ns
    FKS@time[i] <- est$FKS@time
    
    IsolateDetect@dnc[i] <- est$IsolateDetect@nocpt - no.of.cpt
    IsolateDetect@mse[i] <- mean((signal - est$IsolateDetect@fit)^2)
    IsolateDetect@diff <- abs(matrix(est$IsolateDetect@cpt, nrow = no.of.cpt, ncol = length(est$IsolateDetect@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                        ncol = length(est$IsolateDetect@cpt), byr = F))
    IsolateDetect@dh[i] <- max(apply(IsolateDetect@diff, 1, min), apply(IsolateDetect@diff, 2, min))/ns
    IsolateDetect@time[i] <- est$IsolateDetect@time
    
    IsolateDetect_sdll@dnc[i] <- est$IsolateDetect_sdll@nocpt - no.of.cpt
    IsolateDetect_sdll@mse[i] <- mean((signal - est$IsolateDetect_sdll@fit)^2)
    IsolateDetect_sdll@diff <- abs(matrix(est$IsolateDetect_sdll@cpt, nrow = no.of.cpt, ncol = length(est$IsolateDetect_sdll@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                                                      ncol = length(est$IsolateDetect_sdll@cpt), byr = F))
    IsolateDetect_sdll@dh[i] <- max(apply(IsolateDetect_sdll@diff, 1, min), apply(IsolateDetect_sdll@diff, 2, min))/ns
    IsolateDetect_sdll@time[i] <- est$IsolateDetect_sdll@time
    
  }
list(notic = notic, MARS = MARS, t1f = t1f, cpop = cpop, FKS = FKS, 
     IsolateDetect = IsolateDetect, IsolateDetect_sdll = IsolateDetect_sdll)
}

model.wave1 <- list(name = "wave1", cpt.type = "pcwsLinContMean", cpt = c(256, 512, 768, 1024, 1152, 1280, 1344), jump.size = (-1)^(1:7) * 
                      (1:7)/64, n = 1408, start = c(1, 1/256))
lin.SIMR1 = linear_rev.sim.study(model.wave1, 1, seed = 16)


model.wave2 <- list(name = "wave2", cpt.type = "pcwsLinContMean", cpt = (1:99) * 15, jump.size = (-1)^{1:100}, n = 15 * 100, start = c(-1/2, 1/40))
lin.SIMR2 = linear_rev.sim.study(model.wave2, 1, seed = 16)

model.wave3 <- list(name = "wave3", cpt.type = "pcwsLinContMean", cpt = (1:119) * 7, jump.size = (-1)^{1:120}, n = 840, start = c(-1/2, 1/32))
lin.SIMR3 = linear_rev.sim.study(model.wave3, m = 100, sigma = 0.3, seed = 16)

model.wave4 <- list(name = "wave4", cpt.type = "pcwsLinContMean", cpt = (1:9) * 20,
                    jump.size = c(1/6, 1/2,-3/4,-1/3, -2/3,1,1/4,3/4,-5/4), n = 200, start = c(1, 1/32))
lin.SIMR4 = linear_rev.sim.study(model.wave4, m = 100, sigma = 0.3, seed = 16)

model.wave5 <- list(name = "wave5", cpt.type = "pcwsLinContMean", cpt = (1:19) * 50,
                    jump.size = c(-1/16,-5/16,-5/8,1, 5/16,15/32,-5/8,-7/32,-3/4,13/16,5/16,19/32,-1,-5/8,
                                  23/32,1/2,15/16,-25/16,-5/4), n = 1000, start = c(1, 1/32))
lin.SIMR5 = linear_rev.sim.study(model.wave5, m = 100, sigma = 0.6, seed = 16)
