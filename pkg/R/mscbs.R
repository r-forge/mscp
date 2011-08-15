mscbs <-
function(y,f=c(0.01,0.1,1), use.bic=TRUE, MIN.SNPs=3,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,MAX.CHPTS=NA,WCHISQ.CUTOFF=NA,plots=TRUE,CARRIER.RETEST.THRESHOLD=0.01,MIN.REQ.ABSDIFF=NA, MIN.SUFF.ABSDIFF=NA, CHISQ.PVAL.THRESH=0.001){
                                        # For debugging:
                                        # use.bic=TRUE; f=c(0.01,0.1,1); win=1000; Z=NULL; MIN.SNPs=3; ALPHA=0; MAX.CHPTS=NA; WCHISQ.CUTOFF=15; plots=FALSE; CARRIER.RETEST.THRESHOLD=0.01; MIN.REQ.ABSDIFF=NA; MIN.SUFF.ABSDIFF=NA; CHISQ.PVAL.THRESH=0.001

  N=dim(y)[2]
  T=dim(y)[1]
  DELTA = MIN.SNPs/T
  win=T-1
  if(is.na(MAX.CHPTS)) MAX.CHPTS=floor(T/MIN.SNPs)
  if(!use.bic & is.na(WCHISQ.CUTOFF)){
    WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win=T,N,ALPHA)
    cat("MSCBS: weighted chisquare cutoff = ",WCHISQ.CUTOFF,"\n") 
  }
  zlim=c(-1,1)
  
  S = apply(y,2,cumsum)
  y.var<-compute.var(y)
  yhat = matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
  y.r = y-yhat
  SSE.0 = sum(y.r^2)

  if(plots) image.plot(c(1:T),c(1:N),pmax(pmin(yhat,zlim[2]),zlim[1]),zlim=zlim,xlab="position",ylab="sample",main="Fitted");
  
  chpts = c(1,T)
  carriers = rbind(rep(1,N), rep(1,N))
#  bestZ =compute.max.Z.C(y.r,win,y.var,ALPHA,MIN.SNPs)
  bestZ =fcompute.max.Z(y.r,win,y.var,ALPHA,MIN.SNPs)
  best.subchpt= matrix(bestZ$bestchpt,ncol=1,nrow=2)
  best.Z = bestZ$bestZ

  splitnum=0
  chpt.hist = vector("list",MAX.CHPTS)
  bic1.terms = vector("list",MAX.CHPTS)
  bic1 = rep(0,MAX.CHPTS)
  bic2.terms = vector("list",MAX.CHPTS)
  bic2 = rep(0,MAX.CHPTS)
#  yhat.rec = vector("list",MAX.CHPTS)
  SSE = rep(0,MAX.CHPTS)
     
  while(TRUE){    
    if(sum(is.na(best.Z))==length(best.Z)){
        cat("No more good candidates.  MSCBS done.\n",sep="")
        break
    }
    max.Z = max(best.Z,na.rm=TRUE)
    max.region = which.max(best.Z)
    
    if(!use.bic && max.Z<WCHISQ.CUTOFF) {
      cat("Maximum Z-score is ",max.Z,", which does not exceed cutoff of ",WCHISQ.CUTOFF,".  MSCBS done.\n",sep="")
      break
    }

    if(length(chpts)>MAX.CHPTS+2){
      cat("Maximum number of change-points reached.  MSCBS done.\n")
      break
    }
    if(is.na(best.subchpt[1,max.region])){
      cat("Optimal region has no valid change-points.  MSCBS done.\n")      
      break
    }    
    splitnum=splitnum+1
    newchpt = c(best.subchpt[1,max.region],best.subchpt[2,max.region])

    if(use.bic){
        
        if(sum(is.na(newchpt))>0){
            tauhat=c(chpts[1:max.region],newchpt[1], chpts[(max.region+1):length(chpts)])
        } else {
            tauhat=c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])
        }
        
        y.classify = mscbs.classify.bic(y=y,S=S,tauhat=tauhat,y.var=y.var)
        yhat = y.classify$yhat
        y.r = y - yhat
    } else {
        # Classify samples and update yhat.
        y.classify = mscbs.classify(y.r[chpts[max.region]:chpts[(max.region+1)],] , newchpt-chpts[max.region]+1 , y.var , MIN.REQ.ABSDIFF=MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF=MIN.SUFF.ABSDIFF, CHISQ.PVAL.THRESH=CHISQ.PVAL.THRESH)
        y.r.hat = matrix(0,nrow=T, ncol=N)
        y.r.hat[chpts[max.region]:chpts[(max.region+1)],] = y.classify$yhat
        yhat = yhat + y.r.hat
        y.r = y.r - y.r.hat
    }
    

    if(plots) heatmap(pmax(pmin(yhat,zlim[2]),zlim[1]),zlim=zlim, main="Fitted");
        
    if(!is.na(newchpt[2])){  # The added change consistes of two change-points.
      cat("Split ",splitnum,": ",newchpt[1],", ",newchpt[2],", Z-score = ",max.Z,".\n",sep="")
      y.r.L = y.r[chpts[max.region]:(newchpt[1]-1),]
      y.r.M =  y.r[newchpt[1]:(newchpt[2]-1),]
      y.r.R = y.r[newchpt[2]:chpts[max.region+1],]
      
#      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.L = fcompute.max.Z(y.r.L,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.M = compute.max.Z.C(y.r.M,win,y.var,ALPHA,MIN.SNPs)
      bestZ.M = fcompute.max.Z(y.r.M,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = fcompute.max.Z(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      best.Z.new=c(bestZ.L$bestZ, bestZ.M$bestZ, bestZ.R$bestZ)
      best.subchpt.new=cbind(bestZ.L$bestchpt+chpts[max.region]-1, 
        bestZ.M$bestchpt+newchpt[1]-1, 
        bestZ.R$bestchpt+newchpt[2]-1)
    } else { # The added change-point is a singleton.
      newchpt = newchpt[1]
      cat("Split ",splitnum,": ",newchpt,
          ", Z-score = ",max.Z,".\n",sep="")
      y.r.L = y.r[chpts[max.region]:(newchpt-1),]
      y.r.R = y.r[newchpt:chpts[max.region+1],]
      
      
#      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.L = fcompute.max.Z(y.r.L,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = fcompute.max.Z(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      best.Z.new=c(bestZ.L$bestZ, bestZ.R$bestZ)
      best.subchpt.new=cbind(bestZ.L$bestchpt+chpts[max.region]-1, 
        bestZ.R$bestchpt+newchpt-1)
    }
    
    if(max.region>1){ 
      leftpart = best.subchpt[,1:(max.region-1)]
      leftpart.Z=best.Z[1:(max.region-1)]
    }else {
      leftpart = matrix(0,ncol=0, nrow=2)
      leftpart.Z = matrix(0,ncol=0,nrow=0)
    }
    if(max.region+1 <= ncol(best.subchpt)){ 
      rightpart = best.subchpt[,(max.region+1):ncol(best.subchpt)]
      rightpart.Z = best.Z[(max.region+1):length(best.Z)]
    }else  {
      rightpart = matrix(0,ncol=0, nrow=2)
      rightpart.Z=matrix(0,ncol=0,nrow=0)
    }
    
    
    
    if(!use.bic){
        if(length(newchpt)==2){ 
            newcarriers = rbind(y.classify$carriers,y.classify$carriers)
        } else {
            newcarriers = y.classify$carriers
            which.newcarriers = which(newcarriers)
            if(max.region>1) {
                test.left = computeZ.onechange.sample(t(y[chpts[max.region-1]:newchpt,which.newcarriers]),chpts[max.region]-chpts[max.region-1], y.var[which.newcarriers])
                carriers[max.region, which.newcarriers[which(test.left$pval<CARRIER.RETEST.THRESHOLD)]] = 1
            }
            
            if( max.region+2 <= length(chpts)){
                test.right = computeZ.onechange.sample(t(y[newchpt:chpts[max.region+2],which.newcarriers]),chpts[max.region+1]-newchpt, y.var[which.newcarriers])
                carriers[max.region+1, which.newcarriers[which(test.right$pval<CARRIER.RETEST.THRESHOLD)]] = 1
            }
            
        }
        carriers = rbind(carriers[1:max.region,], newcarriers,carriers[(max.region+1):nrow(carriers),])
        phat=NULL
    } else {
        carriers = rbind(rep(1,N), y.classify$carriers,rep(1,N))
        phat = y.classify$phat
    }


    chpt.hist[[splitnum]] = list(chpts=chpts,max.region=max.region,
               newchpt=newchpt,max.Z=max.Z,carriers=carriers)
               
    best.Z = c(leftpart.Z, best.Z.new, rightpart.Z)    
    best.subchpt = cbind(leftpart,best.subchpt.new,rightpart)
    
    chpts = c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])
#    yhat.rec[[splitnum]] = yhat
    
    cat("   Computing BIC: ")
    bic1.terms[[splitnum]] = mbic.ms(S=S,chpts=chpts,sigma=sqrt(y.var),type=1,J=carriers, phat=phat)
    bic1[splitnum] = bic1.terms[[splitnum]]$bic
    bic2.terms[[splitnum]] = mbic.ms(S=S,chpts=chpts,sigma=sqrt(y.var),type=2,J=carriers,phat=phat)
    bic2[splitnum] = bic2.terms[[splitnum]]$bic

    cat(bic1[splitnum]," ",bic2[splitnum],"\n")
    
    SSE[splitnum] = sum(y.r^2)
    cat("   Percent variance explained: ",round(100*(1-(SSE[splitnum]/SSE.0)),digits=2),"\n")
    


 # For debugging: 
    chpts
    rbind(best.Z, best.subchpt)
    
  }

  chpt.hist = chpt.hist[1:splitnum]
  bic2.terms = bic2.terms[1:splitnum]
  bic2 = bic2[1:splitnum]
  bic1.terms=bic1.terms[1:splitnum]
  bic1 = bic1[1:splitnum]
  SSE = SSE[1:splitnum]
   
  if(splitnum>0){ 
    if(use.bic){
        kstar = which.max(bic2)
        if(bic2[kstar]<0){
            chpts.final=rep(0,0)
            yhat.final=matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
            carriers.final = rep(0,0)
            finalsplit=0
        } else {
            chpts.final = c(chpt.hist[[kstar]]$chpts[1:chpt.hist[[kstar]]$max.region],
                            chpt.hist[[kstar]]$newchpt, 
                            chpt.hist[[kstar]]$chpts[(chpt.hist[[kstar]]$max.region+1):length(chpt.hist[[kstar]]$chpts)])  
            chpts.final = chpts.final[2:(length(chpts.final)-1)]
            carriers.final= chpt.hist[[kstar]]$carriers
    #        yhat.final = yhat.rec[[kstar]]     
            yhat.final = mscbs.get.fitted(y=y,tauhat=c(1,chpts.final,nrow(y)),carriers=carriers.final)
            finalsplit=kstar
        }
    } else {
        chpts.final = chpts[2:(length(chpts)-1)]
        finalsplit = splitnum
        carriers.final = carriers
        yhat.final = mscbs.get.fitted(y=y,tauhat=chpts,carriers=carriers)
    }
  } else{
     chpts.final = rep(0,0)
     yhat.final=matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
     finalsplit=0
     carriers.final = rep(0,0)
  }
  
  
  if(plots) heatmap(pmax(pmin(yhat.final,zlim[2]),zlim[1]),zlim=zlim, main="Fitted");
        
  
#  list( chpt.hist=chpt.hist, chpts=chpts.final,
#        yhat=yhat.final, carriers = carriers.final,
#        bic1=bic1, bic1.terms=bic1.terms, bic2=bic2,bic2.terms=bic2.terms, pve = 1-SSE/SSE.0, 
#        finalsplit = finalsplit)

  list( chpt.hist=chpt.hist, chpts=chpts.final,
        yhat=yhat.final, carriers = carriers.final,
        bic=bic2, pve = 1-SSE/SSE.0)


}

