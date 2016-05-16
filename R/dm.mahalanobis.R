dm.mahalanobis <-
function(data,from="mean",p=10,plot=FALSE,v.index=NULL,layout=NULL){
  
  # Initial checks
  if(is.na(match(from,c("mean","median")))){stop('"from" must be mean or median.')}
  if(p>100){stop('"p" must be less than 100(%).')}
  if(trunc(nrow(data)*p/100)==0){stop('"p" is too small (no data to investigate).')}
  if(sum(is.na(apply(data,2,from)))>0){stop('Distance cannot be computed possibly due to nonmetric/missing data.')}
  if(sum(is.na(cov(data)))>0){stop('Covariance cannot be computed possibly due to nonmetric/missing/zero data.')}
  if(!is.null(layout) && (!is.numeric(layout) || length(layout)!=2)){stop('"layout" must be a numeric vector of length 2.')}
  
  # Parameters
  if(isTRUE(plot) && is.null(v.index)){v.index<-1:ncol(data)}
  lo<-function(r){return(c(ceiling(r/5),ceiling(r/ceiling(r/5))))}
  if(isTRUE(plot) && is.null(layout)){layout<-sort(lo(length(v.index)))}
  dp<-function(x){
    if(nchar(abs(x))>1){
      if(trunc(abs(x))==abs(x)){y<-0}
      else{y<-nchar(abs(x))-min(which(abs(x*10^-(1:10))<1))-1}
    }else{y<-0}
    return(ifelse(y>5,2,y))
  }
  digit<-apply(apply(data,c(1,2),dp),2,max)

  # DJL's favorite colors
  col<-c("#FF5A5F","#FFB400","#007A87","#8CE071","#7B0051","#00D1C1","#FFAA91","#B4A76C","#9CA299","#565A5C","#00A04B","#E54C20",
         "thistle","wheat","turquoise","tomato","lightpink","orchid","mediumpurple","lavender","powderblue","gold","lawngreen",
         "mediumspringgreen","palegreen","lemonchiffon","lightsalmon","violet","slateblue","dodgerblue","aquamarine")
  
  # Distance measure
  n<-trunc(nrow(data)*p/100) # n of suspects
  m.dist<-mahalanobis(data,apply(data,2,from),cov(data),tol=1e-25)
  u.order<-order(m.dist,decreasing=TRUE) # dist order in row number
  a.order<-as.numeric(rownames(data)[u.order]) # dist order  in row names
  u.index<-u.order[1:n] # suspect index in row number
  a.index<-as.numeric(rownames(data)[u.index]) # suspect index in row names
  
  # Plot
  if(isTRUE(plot)){
    par(mfrow=layout,mar=c(2,3,4,1)) #c(bottom, left, top, right)
    for(i in v.index){
      boxplot(data[,i],axes=TRUE)
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="whitesmoke",lty=0)
      abline(h=data[u.index,i],col="lightgray")
      boxplot(data[,i],axes=TRUE,main=colnames(data)[i],col=col[i],add=TRUE,cex.lab=1.3)
      text(x=0.65,y=data[u.index,i],labels=paste0("S",a.index," (",format(round(data[u.index,i],digit[i]),nsmall=digit[i]),")"),col="red")
      text(x=1.35,y=fivenum(data[,i]),labels=format(round(fivenum(data[,i]),digit[i]),nsmall=digit[i]))
    }
    par(mfrow=c(1,1))
  }
  
  # Naming
  names(m.dist)<-rownames(data)
  names(u.order)<-c(1:nrow(data))
  names(u.index)<-c(1:n)
  
  # Return results
  results<-list(dist=m.dist,order=u.order,suspects=u.index)
  return(results) 
}
