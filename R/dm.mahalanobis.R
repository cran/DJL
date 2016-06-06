dm.mahalanobis <-
function(data,from="median",p=10,plot=FALSE,v.index=NULL,layout=NULL){
  
  # Effective data
  r.exclude<-which(rowSums(is.na(data))>0); if(length(r.exclude)<1) r.exclude<-NULL
  if(is.null(r.exclude)) data.eff<-data  else data.eff<-data[-r.exclude,]
  
  # Parameter
  n<-nrow(data)
  r.eff<-sum(complete.cases(data))
  
  # Initial checks
  if(is.na(match(from,c("mean","median")))){stop('"from" must be either "mean" or "median".')}
  if(p>100){stop('"p" must be less than 100(%).')}
  if(trunc(nrow(data.eff)*p/100)==0){stop('"p" is too small (no data to be investigated).')}
  if(sum(is.na(apply(data.eff,2,from)))>0){stop('Distance cannot be computed possibly due to nonmetric/missing data.')}
  if(sum(is.na(cov(data.eff)))>0){stop('Covariance cannot be computed possibly due to nonmetric/missing/zero data.')}
  if(!is.null(layout) && (!is.numeric(layout) || length(layout)!=2)){stop('"layout" must be a numeric vector of length 2.')}
  
  # Parameters
  if(plot && is.null(v.index)){v.index<-1:ncol(data.eff)}
  lo<-function(r){return(c(ceiling(r/5),ceiling(r/ceiling(r/5))))}
  if(plot && is.null(layout)){layout<-sort(lo(length(v.index)))}
  dp<-function(x){
    if(nchar(abs(x))>1){
      if(trunc(abs(x))==abs(x)){y<-0}
      else{y<-nchar(abs(x))-min(which(abs(x*10^-(1:10))<1))-1}
    }else{y<-0}
    return(ifelse(y>5,2,y))
  }
  digit<-apply(apply(data.eff,c(1,2),dp),2,max)

  # DJL's favorite colors
  col<-c("#FF5A5F","#FFB400","#007A87","#8CE071","#7B0051","#00D1C1","#FFAA91","#B4A76C","#9CA299","#565A5C","#00A04B","#E54C20",
         "thistle","wheat","turquoise","tomato","lightpink","orchid","mediumpurple","lavender","powderblue","gold","lawngreen",
         "mediumspringgreen","palegreen","lemonchiffon","lightsalmon","violet","slateblue","dodgerblue","aquamarine")
  
  # Distance measure
  n<-trunc(nrow(data.eff)*p/100) # n of suspects
  m.dist<-mahalanobis(data.eff,apply(data.eff,2,from),cov(data.eff),tol=1e-25)
  u.order<-order(m.dist,decreasing=TRUE) # dist order in row number
  a.order<-as.numeric(rownames(data.eff)[u.order]) # dist order in row names
  u.index<-u.order[1:n] # suspect index in row number
  a.index<-as.numeric(rownames(data.eff)[u.index]) # suspect index in row names
  
  # Plot
  if(plot){
    par(mfrow=layout,mar=c(2,3,4,1)) #c(bottom, left, top, right)
    for(i in v.index){
      boxplot(data.eff[,i],axes=TRUE)
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="whitesmoke",lty=0)
      abline(h=data.eff[u.index,i],col="lightgray")
      boxplot(data.eff[,i],axes=TRUE,main=colnames(data.eff)[i],col=col[i],add=TRUE,cex.lab=1.3)
      text(x=0.65,y=data.eff[u.index,i],labels=paste0("S",a.index," (",format(round(data.eff[u.index,i],digit[i]),nsmall=digit[i]),")"),col="red")
      text(x=1.35,y=fivenum(data.eff[,i]),labels=format(round(fivenum(data.eff[,i]),digit[i]),nsmall=digit[i]))
    }
    par(mfrow=c(1,1))
  }
  
  # Returned row numbers
  if(is.null(r.exclude)) r.u.order<-u.order else r.u.order<-c(1:nrow(data))[-r.exclude][u.order]
  if(is.null(r.exclude)) r.u.index<-u.index else r.u.index<-c(1:nrow(data))[-r.exclude][u.index]

  # Naming
  names(m.dist)<-rownames(data.eff)
  names(r.u.order)<-c(1:nrow(data.eff))
  names(r.u.index)<-c(1:n)
  
  # Return results
  results<-list(dist=m.dist,excluded=r.exclude,order=r.u.order,suspect=r.u.index)
  return(results) 
}
