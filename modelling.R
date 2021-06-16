k_selection <- function(model,triadic)
{
  T <- max(which(triadic$event==1))
  k <- 4
  repeat
  {
    answ1 <- StratifiedHazard(model, k)
    answ2 <- StratifiedHazard(model, k+1)
    
    
    l_k <- likelihood(triadic,answ1[[1]],answ1[[2]],answ1[[3]],answ1[[4]],answ1[[5]],answ1[[6]],answ1[[7]],answ1[[8]])
    l_k1 <- likelihood(triadic,answ2[[1]],answ2[[2]],answ2[[3]],answ2[[4]],answ2[[5]],answ2[[6]],answ2[[7]],answ2[[8]])
    
    if(-2*((l_k[[2]] + l_k1[[2]])/2)*(l_k[[1]]-l_k1[[1]]) < 9.49)
    {
      rez <- c(answ1,k)
      break
    }
    if(k==8)
    {
      rez <- c(answ2,9)
      break
    } 
    k <- k+1
  }
  
  return(rez)
}


StratifiedHazard <- function(model, k_value)
{
  max.spon<-max((as.numeric(triadic$stp)-as.numeric(triadic$start))[triadic$strata=="spon"&triadic$event==1])
  max.triad<-max((as.numeric(triadic$stp)-as.numeric(triadic$start))[triadic$strata=="triad"&triadic$event==1])
  max.recip<-max((as.numeric(triadic$stp)-as.numeric(triadic$start))[triadic$strata=="recip"&triadic$event==1])
  max.both <-max((as.numeric(triadic$stp)-as.numeric(triadic$start))[triadic$strata=="both"&triadic$event==1])
  
  bhest_triad <- basehaz(model)
  t <- which(bhest_triad[,3]=="triad")
  baseline_triad <- bhest_triad[t,]
  t <- which(bhest_triad[,3]=="spon")
  baseline_spon <- bhest_triad[t,]
  t <- which(bhest_triad[,3]=="recip")
  baseline_recip <- bhest_triad[t,]
  t <- which(bhest_triad[,3]=="both")
  baseline_both <- bhest_triad[t,]
  
  library(mgcv)
  #triadic closure
  tm.triad<-baseline_triad[,2]
  H.triad<-baseline_triad[,1]
  
  
  triad.approx<-approx(tm.triad,H.triad, n=nrow(baseline_triad))
  tm.triad<-triad.approx$x
  H.triad<-triad.approx$y
  dat<-data.frame(tm.triad,H.triad)
  k <- k_value
  # Show regular spline fit (and save fitted object)
  triad.gam<-gam(H.triad~s(tm.triad,k=k, bs = "cr"),data=dat)
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- smoothCon(s(tm.triad,k=k, bs="cr"),dat,knots=NULL)[[1]]
  F <- mono.con(sm$xp);   # get constraints
  G <- list(X=sm$X,
            C=matrix(0,0,0),
            sp=triad.gam$sp,
            p=sm$xp,
            y=dat$H.triad,
            w=dat$H.triad*0+1)
  G$Ain <- F$A;G$bin <- F$b;G$S <- sm$S;G$off <- 0
  p <- pcls(G);  # fit spline (using s.p. from unconstrained fit)
  H.fit<-Predict.matrix(sm,data.frame(tm.triad=dat$tm.triad))%*%p
  
  newdat.triad<-data.frame(tm.triad=seq(0,max.triad,length=1000))
  triad.pred<-Predict.matrix(sm,data.frame(tm.triad=newdat.triad$tm.triad))%*%p
  
  l=1000
  newdat.triad<-data.frame(tm.triad=seq(0,max.triad,length=l))
  #triad.sd<-as.vector(predict(triad.gam,newdata=data.frame(tm.triad=newdat.triad$tm.triad),type="response", se.fit=TRUE)$se.fit)
  #Variance matrix of parameters
  #GAM model parameters * model matrix
  triad.pred<-Predict.matrix(sm,data.frame(tm.triad=newdat.triad$tm.triad))%*%p
  #Curve variance: X*VP*X^T
  triad.S<-Predict.matrix(sm,data.frame(tm.triad=newdat.triad$tm.triad))%*%triad.gam$Vp%*%t(Predict.matrix(sm,data.frame(tm.triad=newdat.triad$tm.triad)))
  triad.se<-sqrt(diag(triad.S))
  #extract covariance
  triad.cov<-NULL
  for (i in 1:(l-1)){triad.cov[i]<-triad.S[i,i+1]}
  
  
  #reciprocity
  tm.recip<-baseline_recip[,2]
  H.recip<-baseline_recip[,1]
  
  recip.approx<-approx(tm.recip,H.recip, n=nrow(baseline_recip))
  tm.recip<-recip.approx$x
  H.recip<-recip.approx$y
  dat<-data.frame(tm.recip,H.recip)
  
  
  # Show regular spline fit (and save fitted object)
  recip.gam<-gam(H.recip~s(tm.recip,k=k, bs = "cr"),data=dat)
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- smoothCon(s(tm.recip,k=k,bs="cr"),dat,knots=NULL)[[1]]
  F <- mono.con(sm$xp);   # get constraints
  G <- list(X=sm$X,
            C=matrix(0,0,0),
            sp=recip.gam$sp,
            p=sm$xp,
            y=dat$H.recip,
            w=dat$H.recip*0+1)
  G$Ain <- F$A;G$bin <- F$b;G$S <- sm$S;G$off <- 0
  p <- pcls(G);  # fit spline (using s.p. from unconstrained fit)
  H.fit<-Predict.matrix(sm,data.frame(tm.recip=dat$tm.recip))%*%p
  
  newdat.recip<-data.frame(tm.recip=seq(0,max.recip,length=1000))
  recip.pred<-Predict.matrix(sm,data.frame(tm.recip=newdat.recip$tm.recip))%*%p
  
  #l=1000
  newdat.recip<-data.frame(tm.recip=seq(0,max.recip,length=l))
  #recip.sd<-as.vector(predict(recip.gam,newdata=data.frame(tm.recip=newdat.recip$tm.recip),type="response", se.fit=TRUE)$se.fit)
  #Variance matrix of parameters
  #GAM model parameters * model matrix
  recip.pred<-Predict.matrix(sm,data.frame(tm.recip=newdat.recip$tm.recip))%*%p
  #Curve variance: X*VP*X^T
  recip.S<-Predict.matrix(sm,data.frame(tm.recip=newdat.recip$tm.recip))%*%recip.gam$Vp%*%t(Predict.matrix(sm,data.frame(tm.recip=newdat.recip$tm.recip)))
  recip.se<-sqrt(diag(recip.S))
  #extract covariance
  recip.cov<-NULL
  for (i in 1:(l-1)){recip.cov[i]<-recip.S[i,i+1]}
  
  
  #spontaneous
  tm.spon<-baseline_spon[,2]
  H.spon<-baseline_spon[,1]
  
  spon.approx<-approx(tm.spon,H.spon, n=nrow(baseline_spon))
  tm.spon<-spon.approx$x
  H.spon<-spon.approx$y
  dat<-data.frame(tm.spon,H.spon)
  
  
  # Show regular spline fit (and save fitted object)
  spon.gam<-gam(H.spon~s(tm.spon,k=k, bs = "cr"),data=dat)
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- smoothCon(s(tm.spon,k=k,bs="cr"),dat,knots=NULL)[[1]]
  F <- mono.con(sm$xp);   # get constraints
  G <- list(X=sm$X,
            C=matrix(0,0,0),
            sp=spon.gam$sp,
            p=sm$xp,
            y=dat$H.spon,
            w=dat$H.spon*0+1)
  G$Ain <- F$A;G$bin <- F$b;G$S <- sm$S;G$off <- 0
  p <- pcls(G);  # fit spline (using s.p. from unconstrained fit)
  H.fit<-Predict.matrix(sm,data.frame(tm.spon=dat$tm.spon))%*%p
  
  
  newdat.spon<-data.frame(tm.spon=seq(0,max.spon,length=1000))
  spon.pred<-Predict.matrix(sm,data.frame(tm.spon=newdat.spon$tm.spon))%*%p
  
  newdat.spon<-data.frame(tm.spon=seq(0,max.spon,length=l))
  #spon.sd<-as.vector(predict(spon.gam,newdata=data.frame(tm.spon=newdat.spon$tm.spon),type="response", se.fit=TRUE)$se.fit)
  
  #GAM model parameters * model matrix
  spon.pred<-Predict.matrix(sm,data.frame(tm.spon=newdat.spon$tm.spon))%*%p
  #Curve variance: X*VP*X^T
  spon.S<-Predict.matrix(sm,data.frame(tm.spon=newdat.spon$tm.spon))%*%spon.gam$Vp%*%t(Predict.matrix(sm,data.frame(tm.spon=newdat.spon$tm.spon)))
  spon.se<-sqrt(diag(spon.S))
  #extract covariance
  spon.cov<-NULL
  for (i in 1:(l-1)){spon.cov[i]<-spon.S[i,i+1]}
  
  #both
  tm.both<-baseline_both[,2]
  H.both<-baseline_both[,1]
  
  n<-nrow(baseline_both)
  
  
  both.approx<-approx(tm.both,H.both, n=n)
  tm.both<-both.approx$x
  H.both<-both.approx$y

  dat<-data.frame(tm.both,H.both)
  
  
  # Show regular spline fit (and save fitted object)
  both.gam<-gam(H.both~s(tm.both,k=k, bs = "cr"),data=dat)
  # Create Design matrix, constraints etc. for monotonic spline....
  sm <- smoothCon(s(tm.both,k=k,bs="cr"),dat,knots=NULL)[[1]]
  F <- mono.con(sm$xp);   # get constraints
  G <- list(X=sm$X,
            C=matrix(0,0,0),
            sp=both.gam$sp,
            p=sm$xp,
            y=dat$H.both,
            w=dat$H.both*0+1)
  G$Ain <- F$A;G$bin <- F$b;G$S <- sm$S;G$off <- 0
  p <- pcls(G);  # fit spline (using s.p. from unconstrained fit)
  H.fit<-Predict.matrix(sm,data.frame(tm.both=dat$tm.both))%*%p
  
  
  out<-predict(both.gam,newdata=data.frame(tm.both=dat$tm.both),type="response", se.fit=TRUE)

  
  
  l=1000
  newdat.both<-data.frame(tm.both=seq(0,max.both,length=l))
  both.sd<-as.vector(predict(both.gam,newdata=data.frame(tm.both=newdat.both$tm.both),type="response", se.fit=TRUE)$se.fit)
  both.pred<-Predict.matrix(sm,data.frame(tm.both=newdat.both$tm.both))%*%p
  both.S<-Predict.matrix(sm,data.frame(tm.both=newdat.both$tm.both))%*%both.gam$Vp%*%t(Predict.matrix(sm,data.frame(tm.both=newdat.both$tm.both)))
  both.se<-sqrt(diag(both.S))
  both.cov<-NULL
  for (i in 1:(l-1)){both.cov[i]<-both.S[i,i+1]}
  
  l.triad<-diff(triad.pred)/diff(newdat.triad$tm.triad)
  l.spon<-diff(spon.pred)/diff(newdat.spon$tm.spon)
  l.recip<- diff(recip.pred)/diff(newdat.recip$tm.recip)
  l.both<-diff(both.pred)/diff(newdat.both$tm.both)
  s.triad<-sqrt(triad.se[-1]^2+triad.se[-l]^2-2*triad.cov)/diff(newdat.triad$tm.triad)[1]
  s.spon<-sqrt(spon.se[-1]^2+spon.se[-l]^2-2*spon.cov)/diff(newdat.spon$tm.spon)[1]
  s.recip<-sqrt(recip.se[-1]^2+recip.se[-l]^2-2*recip.cov)/diff(newdat.recip$tm.recip)[1]
  s.both<-sqrt(both.se[-1]^2+both.se[-l]^2-2*both.cov)/diff(newdat.both$tm.both)[1]
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  df <- getmode(c(round(summary(triad.gam)$edf),round(summary(spon.gam)$edf),round(summary(recip.gam)$edf)
  ,round(summary(both.gam)$edf)))
  
  return(list(newdat.triad[-1,1],newdat.spon[-1,1],newdat.recip[-1,1],newdat.both[-1,1],
              l.triad[,1],l.spon[,1],l.recip[,1],l.both[,1],
              s.triad, s.spon, s.recip, s.both,df))
}

zero <- function(t1, t_vec, haz)
{
  z1 <- which(abs(t_vec-t1)==min(abs(t_vec-t1)))
  z1 <- haz[z1]
  
  return(ifelse(z1==0,0,1))
}
