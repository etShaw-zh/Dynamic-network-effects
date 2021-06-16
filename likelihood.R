likelihood <- function(triadic,newdat.triad,newdat.spon,newdat.recip,newdat.both,l.triad,l.spon,l.recip,l.both)
{
  
  remove <- NULL
  T <- max(which(triadic$event==1))
  
  t <- which(triadic$strata[1:T] == "spon")
  remove1 <- which(sapply(triadic$stp[t] - triadic$start[t], zero,t_vec=newdat.triad, haz = l.spon )==0)
  remove <- c(remove,t[remove1])
  
  t <- which(triadic$strata[1:T] == "triad")
  remove1 <- which(sapply(triadic$stp[t] - triadic$start[t], zero,t_vec=newdat.spon, haz = l.triad )==0)
  remove <- c(remove,t[remove1])
  
  t <- which(triadic$strata[1:T] == "recip")
  remove1 <- which(sapply(triadic$stp[t] - triadic$start[t], zero,t_vec=newdat.recip, haz = l.recip )==0)
  remove <- c(remove,t[remove1])
  
  t <- which(triadic$strata[1:T] == "both")
  remove1 <- which(sapply(triadic$stp[t] - triadic$start[t], zero,t_vec=newdat.both, haz = l.both )==0)
  remove <- c(remove,t[remove1])
  
  
  seq <- 1:T
  if(length(remove)!=0)
  {
    seq <- seq[-remove]
  }
  pr <- 0
  for(i in seq)
  {
    hazard <- eval(as.name(paste("l.",triadic$strata[i], sep="")))
    time <- eval(as.name(paste("newdat.",triadic$strata[i], sep="")))
    z <- which(abs(time-(triadic$stp[i] - triadic$start[i]))==min(abs(time-(triadic$stp[i] - triadic$start[i]))))
    numerator <- log(hazard[z])
    
    sum <- 0
    
    o <- which(triadic$stp[(i):nrow(triadic)] >= triadic$stp[i] & ((triadic$stp[(i):nrow(triadic)] - triadic$start[(i):nrow(triadic)]) >= (triadic$stp[i] - triadic$start[i])))
    t <- table(triadic$strata[o+(i-1)])
    
    
    hazard1 <- l.spon
    time1 <- newdat.spon
    n1 <- which(names(t)=="spon")
    if(length(n1)!=0)
    {
      temp <- which(abs(time1-(triadic$stp[i] - triadic$start[i]))==min(abs(time1-(triadic$stp[i] - triadic$start[i]))))
      sum <- sum + hazard1[temp]*t[n1]
    }
    
    
    hazard1 <- l.triad
    time1 <- newdat.triad
    n1 <- which(names(t)=="triad")
    if(length(n1)!=0)
    {
      temp <- which(abs(time1-(triadic$stp[i] - triadic$start[i]))==min(abs(time1-(triadic$stp[i] - triadic$start[i]))))
      sum <- sum + hazard1[temp]*t[n1]
    }
    
    hazard1 <- l.recip
    time1 <- newdat.recip
    n1 <- which(names(t)=="recip")
    if(length(n1)!=0)
    {
      temp <- which(abs(time1-(triadic$stp[i] - triadic$start[i]))==min(abs(time1-(triadic$stp[i] - triadic$start[i]))))
      sum <- sum + hazard1[temp]*t[n1]
    }
    
    hazard1 <- l.both
    time1 <- newdat.both
    n1 <- which(names(t)=="both")
    if(length(n1)!=0)
    {
      temp <- which(abs(time1-(triadic$stp[i] - triadic$start[i]))==min(abs(time1-(triadic$stp[i] - triadic$start[i]))))
      sum <- sum + hazard1[temp]*t[n1]
    }
    
    pr <- pr + (numerator-log(sum))
  }
  return(list(pr/length(seq),length(seq)))
}
