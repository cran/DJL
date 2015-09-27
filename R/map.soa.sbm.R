map.soa.sbm <-
function(xdata,ydata,date,rts,orientation,sg="ssm"){
  # Subset index
  till<-function(x,y){
    t<-0
    while(x[t+1]<=y&&t<nrow(x)){t<-t+1}
    return(t)
  }
  
  # Parameters
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  
  # Sort data ascending order
  x<-matrix(c(xdata[order(date),]),ncol=m)
  y<-matrix(c(ydata[order(date),]),ncol=s)
  d<-matrix(c(date[order(date),]),ncol=1)
  
  # max map size
  c<-nrow(unique(d)) 
  ud<-unique(d)
  
  # map frame
  fanta<-matrix(c(NA),nrow=n,ncol=c);colnames(fanta)<-ud
  
  # generate the map
  for(i in 1:c){
    # subset data
    e<-till(d,ud[i])
    x_s<-matrix(x[1:e,],nrow=e)
    y_s<-matrix(y[1:e,],nrow=e)
    
    # run distance measure
    dj<-dm.sbm(x_s,y_s,rts,orientation="n",se=0,sg="ssm")
    
    # soa set
    soa<-which(round(dj$eff,8)==1)
    #soa<-intersect(which(round(dj$eff,8)==1), which(cbind(dj$xslack,dj$yslack)==0))
    
    # fill the map
    j<-sum(soa>0)
    q<-1
    for(k in 1:j){
      if(ud[i]==ud[1]){fanta[k,1]<-soa[k]}
      else{
        l<-which(fanta[,i-1]==soa[k])
        if(length(l)>0){fanta[l,i]<-soa[k]}
        else{
          p<-n
          while(is.na(fanta[p,i-1])){p<-p-1}
          fanta[p+q,i]<-soa[k]
          q<-q+1
        }
      }
    }
  }
  fanta<-fanta[1:(p+q-1),]
  rownames(fanta)<-na.omit(unique(c(fanta)))
  print(fanta)
}
