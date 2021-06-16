#Upload required function
source("data preparation.R")
source("strata.R")
source("likelihood.R")
source("modelling.R")
#Data simulation
#number of nodes
n=20
#number of events
E=1000
#model parameters
lambda.spon = 0.2
lambda.triad = 0.5
lambda.recip = 0.8
lambda.both = 1.2
info <- cbind(rep(1:n,each=n),rep(1:n,n),0,0)
t <- which(info[,1] == info[,2])
info <- info[-t,]
info <- as.data.frame(info)
info <- cbind(info,rep("spon",nrow(info)))
info <- cbind(info,rep(lambda.spon,nrow(info)))
colnames(info) <- c("from","to","triad","recip","state","lambda")
info[,5] <- as.character(info[,5])
simdat<- matrix(0,E,4)
simdat <- as.data.frame(simdat)
simdat <- unname(simdat)
st=0
for (i in 1:E)
{
  parameters <- info[,6]
  #time until event
  tm<-rexp(1,sum(parameters))
  #which link
  link <- rmultinom(1,1, prob = c(parameters)/sum(parameters))
  id <- which(link==1)
  #id <- kk[i]
  st<-st+tm
  hlp<-c(info[id,1:2],st,info[id,6])
  hlp
  simdat[i,] <-hlp
  s <- info[id,1]
  r <- info[id,2]
  #check if it creates reciprocity
  if(info[info[,1]==s & info[,2]==r,4]==0 & info[info[,1]==r & info[,2]==s,4]==0)
  {
    info[info[,1]==r & info[,2]==s,4] <- 1
    #check if it is triadic
    if(info[info[,1]==r & info[,2]==s,3] > 0)
    {
      info[info[,1]==r & info[,2]==s,6] <- lambda.both
    } else{
      info[info[,1]==r & info[,2]==s,6] <- lambda.recip
    }
  }
  #check if it creates an open triad
  #if it happened before:
  if(i > 2)
  {
    prev <- which(unlist(simdat[1:(i-1),1])==s & unlist(simdat[1:(i-1),2])==r)
    prev <- ifelse(length(prev)>0,prev[length(prev)]+1,1)
    #find from whom sender receives emails
    t <- which(unlist(simdat[prev:i,2])==s)
    t <- prev + t - 1
    l <- unlist(simdat[t,1])
    l <- unique(l)
    if(length(l)!=0)
    {
      for(j in 1:length(l))
      {
        if(l[j]!=r)
        {
          #if it was non triadic
          if(info[info[,1]==r & info[,2]==l[j],3] == 0)
          {
            #if it was non reciprocal
            if(info[info[,1]==r & info[,2]==l[j],4] == 0)
            {
              info[info[,1]==r & info[,2]==l[j],6] <- lambda.triad
            } else{
              info[info[,1]==r & info[,2]==l[j],6] <- lambda.both
            }
          }
          info[info[,1]==r & info[,2]==l[j],3] <- info[info[,1]==r & info[,2]==l[j],3] + 1
        }
      }
    }
  }
  #-------------Define link status------------
  #check if link was reciprocal
  if(info[info[,1]==s & info[,2]==r,4]==1)
  {
    info[info[,1]==s & info[,2]==r,4] <- 0
    #check if it was triadic
    if(info[info[,1]==s & info[,2]==r,3] > 0)
    {
      info[info[,1]==s & info[,2]==r,3] <- 0
    }
    info[info[,1]==s & info[,2]==r,6] <- lambda.spon
  } else{
    #check if it closes a triad
    if(info[info[,1]==s & info[,2]==r,3] > 0)
    {
      info[info[,1]==s & info[,2]==r,3] <- 0
      info[info[,1]==s & info[,2]==r,6] <- lambda.spon
    }
  }
}
simdat <- cbind(simdat)
#data preparation
#include censored events
data <- e2(simdat[,c(3,1,2)],n)
#assign events to corresponding strata
triadic <- triad_f(data)

#---------------Modelling--------------------
library(survival)
surv <- Surv(as.numeric(triadic$stp)-as.numeric(triadic$start), as.numeric(triadic$event))
model <- coxph(surv ~ 1 + strata(triadic$strata),robust=TRUE)
answ <- k_selection(model,triadic)

plot(answ[[1]],answ[[5]],type="l",col="green",ylim=c(0,2), bty="n", xlab='Time', ylab='Baseline hazard',xaxt="n",yaxt='n',cex.lab=1.7,lwd=2.5)
lines(answ[[2]],answ[[6]],type="l",col="red",lwd=2.5)
lines(answ[[2]],as.vector(answ[[6]])+1.96*as.vector(answ[[10]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[3]],as.vector(answ[[7]])-1.96*as.vector(answ[[11]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[3]],answ[[7]],type="l",col="blue",lwd=2.5)
lines(answ[[4]],answ[[8]],type="l",col="yellow",lwd=2.5)
lines(answ[[4]],as.vector(answ[[8]])-1.96*as.vector(answ[[12]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[4]],as.vector(answ[[8]])+1.96*as.vector(answ[[12]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[2]],as.vector(answ[[6]])-1.96*as.vector(answ[[10]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[3]],as.vector(answ[[7]])+1.96*as.vector(answ[[11]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[1]],as.vector(answ[[5]])-1.96*as.vector(answ[[9]]),type="l",col="black",lwd=2.5,lty=2)
lines(answ[[1]],as.vector(answ[[5]])+1.96*as.vector(answ[[9]]),type="l",col="black",lwd=2.5,lty=2)
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)
legend("topright", legend=c("Spontaneous", "Triadic effect","Reciprocity","Both"),
col=c("red", "green","blue","yellow"), lty=1, cex=1.5,bty = "n", lwd=2.5)

