e2 <- function(e,n)
{
  n.events <- nrow(e)
  #time moments from 0 to event time
  sender <- e[,2]
  receiver <- e[,3]
  start<-rep(0,n.events)
  stp<-e[1:n.events,1]
  #event happened - 1, not happened - 0
  event<-rep(1,n.events)
  for (i in 1:(n.events - 1)){
    #checks if the same event happens more than once
    same.link<-(e[(i+1):n.events,2]==e[i,2])&(e[(i+1):n.events,3]==e[i,3])
    #if is true (i.e., happens more than once)
    if (sum(same.link)>0){
      #index of this link
      ind<-((i+1):n.events)[same.link][1]
      #start of link ind
      start[ind]<-e[i,1]
    } else{
      #it happened, but only once, it means, that it is still "active"
      start<-c(start,e[i,1])
      #observing until this time point
      stp<-c(stp,(max(e[,1])+min(diff(e[,1]))/2))
      #and it not active
      event<-c(event,0)
      sender <- c(sender,e[i,2])
      receiver <- c(receiver,e[i,3])
    }
  }
  #all possible pairs
  for (s in 1:n){
    for (r in 1:n){
      if (s!=r){
        #if this link is never present
        if (sum((e[1:n.events,2]==s) & (e[1:n.events,3]==r))==0){
          #its beginning is here
          start<-c(start,0)
          #observation end
          stp<-c(stp,(max(e[,1])+min(diff(e[,1]))/2))
          #it is not active 
          event<-c(event,0)
          sender <- c(sender,s)
          receiver <- c(receiver,r)
        }
      }
    }
  }
  data <- cbind(sender,receiver,start,stp,event)
  data <- as.data.frame(data)
  return(data)
}