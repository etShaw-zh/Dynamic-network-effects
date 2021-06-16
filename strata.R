triad_f <- function(data)
{
  #number of actors
  n <- length(unique(c(data[,1],data[,2])))
  #number of events
  n.events <- nrow(data)
  
  triadic <- as.data.frame(data)
  triadic <- cbind(triadic, rep("spon",n.events))
  triadic[,6] <- as.character(triadic[,6])
  colnames(triadic)[6] <- "strata"
  
  info <- cbind(rep(1:n,each=n),rep(1:n,n),0,0,0,0)
  t <- which(info[,1] == info[,2])
  info <- info[-t,]
  info <- as.data.frame(info)
  info <- cbind(info,rep("spon",nrow(info)))
  colnames(info) <- c("from","to","triad","time","closed_triad","recip","state")
  info[,7] <- as.character(info[,7])
  for(i in 1:nrow(data))
  {
    s <- data$sender[i]
    r <- data$receiver[i]
    #if event happens
    if(data$event[i]==1)
    {
      #check if it creates reciprocity
      if(info[info[,1]==s & info[,2]==r,6]==0 & info[info[,1]==r & info[,2]==s,6]==0)
      {
        info[info[,1]==r & info[,2]==s,6] <- 1
        #check if it is triadic
        if(info[info[,1]==r & info[,2]==s,3] > 0)
        {
          new <- c(r,s,info[info[,1]==r & info[,2]==s,4],data[i,4],0,"triad")
          triadic <- rbind(triadic,new)
        } else{
          new <- c(r,s,info[info[,1]==r & info[,2]==s,4],data[i,4],0,"spon")
          triadic <- rbind(triadic,new)
        }
        info[info[,1]==r & info[,2]==s,4] <- data[i,4]
      }
      
      #check if it creates an open triad
      #if it happened before:
      if(i > 2)
      {
        prev <- which(data[1:(i-1),1]==s & data[1:(i-1),2]==r)
        prev <- ifelse(length(prev)>0,prev[length(prev)]+1,1)
        #find from whom sender receives emails
        t <- which(data[prev:i,2]==s)
        t <- prev + t - 1
        l <- data[t,1]
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
                if(info[info[,1]==r & info[,2]==l[j],6] == 0)
                {
                  new <- c(r,l[j],info[info[,1]==r & info[,2]==l[j],4],data[i,4],0,"spon")
                  triadic <- rbind(triadic,new)
                } else{
                  new <- c(r,l[j],info[info[,1]==r & info[,2]==l[j],4],data[i,4],0,"recip")
                  triadic <- rbind(triadic,new)
                }
                info[info[,1]==r & info[,2]==l[j],4] <- data[i,4]
              } 
              info[info[,1]==r & info[,2]==l[j],3] <- info[info[,1]==r & info[,2]==l[j],3] + 1 
              
            }
          }
        }
      }
      
      
    }
    
    #-------------Define link status------------
    #check if link was reciprocal
    if(info[info[,1]==s & info[,2]==r,6]==1)
    {
      #check if it was triadic
      if(info[info[,1]==s & info[,2]==r,3] > 0)
      {
        info[info[,1]==s & info[,2]==r,3] <- 0
        triadic$strata[i] <- "both"
      } else{
        triadic$strata[i] <- "recip"
        
      }
      info[info[,1]==s & info[,2]==r,6] <- 0 
      
    } else{
      #check if it closes a triad
      if(info[info[,1]==s & info[,2]==r,3] > 0)
      {
        info[info[,1]==s & info[,2]==r,3] <- 0
        triadic[i,3] <- info[info[,1]==s & info[,2]==r,4]
        triadic$strata[i] <- "triad"
      }
    }
    triadic[i,3] <- info[info[,1]==s & info[,2]==r,4]
    info[info[,1]==s & info[,2]==r,4] <- data[i,4] 
    
  }
  colnames(triadic)[1:2] <- c("sender","receiver")
  triadic$start <- as.numeric(triadic$start)
  triadic$stp <- as.numeric(triadic$stp)
  return(triadic)
}
