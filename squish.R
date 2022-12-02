### Take a matrix and conversion list and squish them (newer function)
squish<-function(inexp,convlist,verbose=FALSE,method="sum"){
  inexp<-inexp[!is.na(convlist),]
  egs<-convlist[!is.na(convlist)]
  
  outexp<-matrix(0,nrow=length(unique(egs)),ncol=ncol(inexp))
  rownames(outexp)<-unique(egs)
  colnames(outexp)<-colnames(inexp)
  
  
  # Sum values of Rows mapping the same Gene
  if(method=="sum"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-which(egs==eghere)
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,sum,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  # Average values of Rows mapping the same Gene
  if(method=="average"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,mean,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
    # Average values of Rows mapping the same Gene using median
  if(method=="median"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,median,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
 # Keep only max value of Rows mapping the same Gene (most expressed probe)
  if(method=="highest"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,max,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  #TODO: keep the highest variance, keep the most expressed, etc.
  
  
  
  return(outexp)
}
