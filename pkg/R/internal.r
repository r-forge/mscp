

BIC.means.nop  <-function(bX,n=NULL){
    m=length(bX)
    if(length(n)==1){
        n = rep(n,m)
    }

    term1 = n*bX^2/2 - log(n)/2
    ord = order(term1, decreasing=TRUE)

    J.size = c(1:m)
    c = cumsum(bX[ord]^2)/J.size
    
    bic.term1 = cumsum(term1[ord])-0.5*J.size-0.5*log(c)
    bic = bic.term1 - log(J.size)/2
    
    bic = c(0,bic)
    calls = which.max(bic)-1
    
    if(calls==0){
        Jhat = rep(0,0)
    } else {
        Jhat = ord[1:calls]
    }   
    
    list(model.size = c(0:m), bic=bic, calls = calls, Jhat =Jhat )
    

}


BIC.means <-function(bX,n=NULL,p=NULL,plots=TRUE, plot.heatmap=FALSE, plot.mfrow=NULL, truth=NULL){
    m=length(bX)
    if(length(n)==1){
        n = rep(n,m)
    }

    term0 = n*bX^2/2
    term1 = n*bX^2/2 - log(n)/2
    ord = order(term1, decreasing=TRUE)

    if(is.null(p)){
        p.min = 1/(2*m)
        p.max = 1-1/(2*m)
        p = seq(p.min, p.max, 1/(2*m))
    }
    
    
    pmat = matrix(nrow=m, ncol=length(p),data=p,byrow=TRUE)
    term1mat = matrix(nrow=m, ncol=length(p),data=term1[ord], byrow=FALSE)
    term2 = pmax(term1mat - log((1-pmat)/pmat),0)
    term3 = log(1-pmat)
    tbic = apply(term2+term3,2,sum) 


    J.size = c(1:m)
    c = cumsum(bX[ord]^2)/J.size
    
    bic.term1 = cumsum(term1[ord])-0.5*J.size-0.5*log(c)
    
    J.vec = matrix(data=J.size,nrow=m,ncol=1)
    p.vec = matrix(data=p,nrow=1,ncol=length(p))
    bic.term2 = J.vec%*%log(p.vec)+(m-J.vec)%*%log(1-p.vec) 
    
    bic = matrix(data=bic.term1,nrow=m,ncol=length(p),byrow=FALSE) + bic.term2 
            - matrix(data=log(J.size)/2, nrow=m,ncol=length(p),byrow=FALSE)
    bic.max = apply(bic,2,max)
    
    
    p = c(0,p)
    tbic= c(0,tbic)
    tbic.p.ind = which.max(tbic)
    if(tbic.p.ind>1) {
        tbic.calls = sum(term2[,tbic.p.ind-1]>0)
    } else {
        tbic.calls = 0
    }
    
    bic=cbind(rep(0,m),bic)
    bic.max = c(0,bic.max)
    bic.p.ind =  which.max(bic.max)
    if(bic.p.ind>1) {
        bic.calls = which.max(bic[,bic.p.ind])
    } else {
        bic.calls = 0
    }
    
    
    bic.nop = BIC.means.nop(bX,n)
    
    if(plots){
        if(is.null(plot.mfrow)){
            if(plot.heatmap){
                par(mfrow=c(3,1))
            }else{
                par(mfrow=c(2,1))
            }
        } else {
            par(mfrow=c(plot.mfrow,1))
        }
        plot(sqrt(n)*abs(bX[ord]), type="b", ylim=c(0,max(sqrt(n)*abs(bX[ord]))), xlab="Sorted rank", main="Z-scores, ranked by absolute value",ylab="z-score", cex.main=1.5, cex.lab=1.5)
#        plot(term1[ord], type="l", ylim=c(min(term1),max(term0)), xlab="Sorted rank", ylab="sqrt(n)*(absolute mean)", cex.lab=1.5)
#        points(term0[ord],col="black")
        segments(tbic.calls,0, tbic.calls,max(term0), col="blue",lty=2)
        segments(bic.calls,0, bic.calls,max(term0), col="red",lty=3)
        segments(bic.nop$calls,0,bic.nop$calls,max(term0), col="orange",lty=4)
        legend(x="topright",lty=c(2,3,4),col=c("red","blue","orange"),legend=c("BIC","simplified BIC","BIC w/ flat J prior"))
        
        
        
        plot(p,tbic, type="l", col="blue", ylim=c(min(bic.max),max(tbic)), main="BIC and simplified BIC",cex.main=1.5, xlab="p", ylab="BIC", cex.lab=1.5)
        lines(p,bic.max, col="red")
        segments(p[tbic.p.ind],min(bic.max),p[tbic.p.ind],max(tbic),col="blue",lty=2)
        segments(p[bic.p.ind],min(bic.max),p[bic.p.ind],max(tbic),col="red",lty=3)
        if(!is.null(truth)) segments(truth,min(bic.max),truth,max(tbic),col="black",lty=1)
        abline(0,0)
        legend(x="bottom",lty=c(1,1,1),col=c("red","blue"),legend=c("BIC","simplified BIC"))
        
        if(plot.heatmap){
    #        cols=rainbow(m)    
    #        plot(p, bic[1,], type="l", col=cols[i],ylim=c(min(bic),0))
    #        for(i in 2:m){
    #            lines(p,bic[i,], col=cols[i])
    #        }
    #        lines(p,bic.max)
            
            image.plot(p,J.size,t(bic), xlab="p", ylab="|J|", cex.lab=1.5, main="Heatmap of BIC", cex.main=1.5 )
        }

    }

    if(bic.calls>0){
        bic.J = ord[1:bic.calls]
    } else {
        bic.J = rep(0,0)
    }
    
    
    if(tbic.calls>0){ 
        tbic.J = ord[1:tbic.calls]
    } else {
        tbic.J = rep(0,0)
    }


    list(p=p,bic=bic,bic.max = bic.max, tbic.max=tbic,
        bic.calls = bic.calls, tbic.calls=tbic.calls, 
        bic.phat = p[bic.p.ind], tbic.phat = p[tbic.p.ind], 
        bic.J=bic.J, tbic.J=tbic.J)
}

Change.Points.R <- function(ix, T, lookup){

    from.row <- (ix[1:lookup] %% T)
    from.col <- as.integer( ix[1:lookup] / T) + 1
    
    v1 <- from.row + 1
    v2 <- from.col + from.row
    return(cbind(v1, v2))
}

Chisq.Contrib.R <- function(y, T, chpts){

# chpts = [st, ed], changed segment runs from st+1 to ed.

    dfnum <- 1;
    dfden <- T-2;

    N <- nrow(y)
    ret.m <- matrix(nrow=nrow(chpts), ncol=N, data=0)
    # chisq holds the contribution of each sample to the aberration.

    n.segs <- nrow(chpts)

    for(i in 1:n.segs){
        st <- chpts[i,1];
        ed <- chpts[i,2];
        w <- ed-st; # Changed 11/4, previously w<-ed-st+1.

        for(j in 1:N){
            mn1 <- sum(y[j,])/T
            SST <- sum( (y[j,]- mn1)^2 );
            SSb <- (( sum(y[j,(st+1):ed])-w*mn1)^2) / (w*(1-w/T));
            SSw = SST-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/dfden);
        }

    }
    return(ret.m)
}

Chisq.Contrib.C <- function(y, T, chpts){
    res1 <- .Call("ChisqContrib", as.vector(y), T, dim(y), as.integer(chpts[,1]), as.integer(chpts[,2])  )
    return(matrix(data=res1, nrow=nrow(chpts))  )
}

Chisq.Contrib.test <- function(){
    y <- matrix(runif(100*20), ncol=20)
    T <- 5

    chpts <- matrix(runif(2*20), ncol=2)
    chpts[,2] <- rowSums(chpts)
    chpts[,1] <- as.integer(chpts[,1] * 15) + 1
    chpts[,2] <- as.integer(chpts[,2] * 15) + 3

    chpts[chpts > ncol(y)] <- ncol(y)
    
    
    res1 <- Chisq.Contrib.R(y, T, chpts)
    res2 <- Chisq.Contrib.C(y, T, chpts)

    print( range(abs(res1-res2))  )

    par(mfrow=c(1,3))
    image(res1)
    image(res2)
    image(res1-res2)
}


Remove.Overlap.R <- function(chpts){

    lookup <- nrow(chpts)
    
    ret.v <- rep(FALSE, lookup)
    ret.v[1] <- TRUE
    for(i in 2:lookup){
        overlap <- 0;
        st <- chpts[i,1];
        ed <- chpts[i,2];
        for(j in which(ret.v)){
            # check for overlap.
            if(     (st<=chpts[j,1] && ed>=chpts[j,2]) || 
                (st>=chpts[j,1] && st<=chpts[j,2]) || 
                (ed>=chpts[j,1] && ed<=chpts[j,2]) )
            {
                overlap <- 1;
                break;
            }
        }
        if(overlap == 0){
            ret.v[i] <- TRUE
        }
            
    }

    return(ret.v)
}


Remove.Overlap.C <- function(chpts){
#This is a C version of the function above.
    res1 <- .Call("RemoveOverlap", chpts[,1], chpts[,2]);
    
    i2 <- rep(FALSE, nrow(chpts))
    i2[res1 == 1] <- TRUE
    return(i2)
}

Remove.Overlap.test <- function(){
    par(mfrow=c(1,3))

    y <- matrix(runif(2*20), ncol=2)
    y[,2] <- rowSums(y)

    i2a <- Remove.Overlap.R(y)
    i2b <- Remove.Overlap.C(y)

    

    print( sum(i2a != i2b) )

    i2a <- matrix(as.integer(i2a), nrow=1)
    i2b <- matrix(as.integer(i2a), nrow=1)

    image(i2a)
    image(i2b)
    image(i2a - i2b)
    
}




ComputeZ.R <- function(Y, T, win, ALPHA){

#  This function does not take in psidot0 as a parameter, so 
#  after calling Z should be modified to Z = Z - N * psidot0
#  and then the final operation  'Z=Z/sqrt(psidotdot0*N);'
#  should be performed.

    dfnum <- 1;
    dfden <- T-2;



    log.alpha <- log(ALPHA)
    
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    

    N <- dim(y)[1]

    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;


    for(i in 1:N){

        
        S=cumsum(y[i,]);
        SST = sum((y[i,]-S[T]/T)^2);
        
        for(k in 1:win){
            diff1 <- S[(k+1):T]-S[1:(T-k)]
                SSb = k*(diff1/k-S[T]/T)^2;
                SSb = SSb + (T-k)*((S[T]-diff1)/(T-k)-S[T]/T)^2;
                SSw = SST-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
        }
        
        Z <- Z+g(U);  # Z[t,k]: change starting at t+1, ending at t+k.
    }
    return(Z)
}

computeZ.onechange<-function(y,t,y.var){
  T =ncol(y)
# cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y),"\n")

  St = apply(y[,1:t,drop=FALSE],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  sumZ=sum(Zsq)
  pval = 1-pchisq(sumZ,nrow(y))
#  cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y)," sumZ=", sumZ,", pval=", pval,"\n")
  # Looks like the Z-values are immense, and rarely ever accepts the null.
#  cat(Zsq)
#  cat("\n\n")
  
  list(sumZ=sumZ, pval=pval)
}

computeZ.onechange.sample<-function(this.y,t,this.y.var){
# unlike computeZ.onechange, this function finds the carriers.

  T =ncol(this.y)
  N = nrow(this.y)
 # cat("here!! T=", T,", t=", t,", nrow(this.y)=", nrow(this.y),"\n\n")


  St = apply(matrix(data=this.y[,1:t],nrow=N,ncol=t),1,sum)
  ST = apply(matrix(data=this.y, nrow=N,ncol=T),1,sum)
  Zsq = (St - (t/T)*ST)^2/(this.y.var*t*(1-t/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq, pval=pval)
}

computeZ.squarewave.sample<-function(this.y,seg,y.var){
# assumes that changed region is seg[1]+1 to seg[2]
  T =ncol(this.y)
  Ss = apply(as.matrix(this.y[,1:seg[1]],nrow=nrow(this.y),ncol=seg[1]),1,sum)
  St = apply(this.y[,1:seg[2]],1,sum)
  ST = apply(this.y,1,sum)
  k=seg[2]-seg[1]
  Zsq = (St-Ss - (k/T)*ST)^2/(y.var*k*(1-k/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq=Zsq, pval=pval)
}


ComputeZ.C <- function(y, T, win, ALPHA){
    res1 <- .Call("ComputeZ", as.vector(y), T, win, dim(y), ALPHA);

    Z <- matrix(nrow=T, ncol=win, data=res1)
    return(Z)
}


computeZ.test <- function(){
    par(mfrow=c(1,3))


    y <- matrix(rnorm(100*20), nrow=20)

    T = 100;
    win = 20;

 #   system.time(z <- ComputeZ.R(y, T, win, 0))
    system.time(z2 <- ComputeZ.C(y, T, win, 0))
    
    print(range(abs(z-z2))  )

    image(z)
    image(z2)
    image(z-z2)
    
}

flatten <- function(z){
# takes the values of matrix z one COLUMN at a time
# and put them in a vector.
  n=prod(dim(z))
  z.flat = matrix(data=z,nrow=n, ncol=1)
}


plotChpts <- function(msscan.retval){

    plot(   c(0, ncol(msscan.retval$yhat)),
        c(0, nrow(msscan.retval$yhat)),
        xlab="sample index",
        ylab="SNP index",
        col="white",
        main="")
    
    for(j in 1:nrow(msscan.retval$yhat)){
        v1 <- msscan.retval$yhat[j,-1]
        v2 <- msscan.retval$yhat[j,-ncol(msscan.retval$yhat)]
        
        i2 <- which(v2 != v1)

        if(length(i2) == 0) next

        for(k in i2){
            lines(c(k, k+1, k+1, k,k), c(j,j,j+1,j+1,j), lty="solid")
        }
    }
        
}

getCutoffMultisampleWeightedChisq <- function(pval,m,delta,win,N,alpha){
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES = 0.1*pval
    currb = 1
    prevsmallerb = currb
    prevlargerb = 200
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
#        cat("pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
#        cat("currb = ",currb,"\n")
        currpval = pvalueMultisampleWeightedChisq(currb,m,delta,win,alpha,N);
    }    
    currb
}


pvalueMultisampleWeightedChisq<-function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b*sqrt(beta)/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
}



pvalueMultisampleWeightedChisqold<-function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
}


mllr<-function(u,r){
    u^2/2 + log((1-r)*exp(-u^2/2) + r)
}
mllrdot<-function(u,p0){
    p0*u/((1-p0)*exp(-u^2/2) + p0)
}

wchi<-function(u,p0){
    alpha=(1-p0)/p0
    eterm = exp(-u^2/2)*alpha
    u^2/(eterm+1)
}

wchidot<-function(u,p0){
    alpha=(1-p0)/p0
    eterm = exp(-u^2/2)*alpha
    (2*u*(eterm+1)+eterm*u^3/2)/(1+eterm)^2
}



wschi<-function(X,r){
    expterm=exp(-X^2/2)
    alpha=(1-r)/r
    stat = sum(X^2/(alpha*expterm+1))
    stat
}




nuFunction<-function(x){
    ((2/x)*(pnorm(x/2)-.5))/((x/2)*pnorm(x/2)+dnorm(x/2))
}

psiNormal<-function(theta,g,p0){

    INTLIM.THRESH=0.001
    psidotbot.int =  function(u){ exp(theta*g(u,p0)-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
   
    log(psidotbot)
}


computeLocalIncrementMean<-function(g,gdot,theta,p0, LOWER=-50, UPPER=50){
    
    psi = psiNormal(theta,g,p0)
    
    integrand<-function(z){
        gdot(z,p0)^2*exp(theta*g(z,p0)-psi-.5*z^2)
    }
    
#    LOWER=-50; UPPER=50
#    temp=seq(LOWER,UPPER,.1)
#    plot(temp, integrand(temp))
    
    
    integral = integrate(integrand, lower=LOWER,upper=UPPER)
    
    .5*theta^2*integral$value/sqrt(2*pi)
    
}



getCutoffMultisample <- function(pval,m,m0,m1,N,g,gdot,p0){
    cat("\nComputing p-value cut off for m=",m,", m0=",m0,", m1=",m1,", N=",N,", p0=",p0,", pvalue=",pval,".\n",sep="")

    func<-function(z){ g(z,p0)*exp(-z^2/2)}
    Eg = integrate(func,lower=-10, upper=10)
    Eg = Eg$value/sqrt(2*pi)
    func<-function(z){ (g(z,p0)-Eg)^2*exp(-z^2/2)}
    Vg = integrate(func, lower=-10, upper=10)
    Vg = Vg$value/sqrt(2*pi)

    THRES = 0.1*pval
    currb = N*Eg+1*sqrt(Vg*N)
    prevsmallerb = currb
    prevlargerb = N*Eg+100*sqrt(Vg*N)
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
        cat("     getCutoffMultisample: pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
        cat("     currb = ",currb,"\n")
        currpval = pvalueMultisample(currb,m,m0,m1,N,g,gdot,p0)
    }    
    currb
}


pvalueMultisample<-function(b.3,m,m0,m1,N,g,gdot,p0){
# b: threshold
# m: total length
# m0: smallest window size
# m1: largest window size
# N: number of sequences
# g: the transformation function
# p0: weight factor.
# b/N must be larger than psidot(theta=0,b=0), which is the mean under the null.
# 
# Note: pvalueMultisample(b.3,..mllr,mllrdot,p0=1) should be equal to pvalueSumChisqNew(2*b.3,...)
# 
# Added 8/15: pvalueMultisample(b.3,...,g,gdot,p0) is the probability P(max_{s,t} \sum_{i=1}^N g(U_{ist}^2,p0) > b.3)
# 
    
    theta0=0.1
    gr <-function(u){ g(u,p0)}
    numtries=0
    MAX.TRIES=100
    while(numtries<MAX.TRIES){
        theta0 = sqrt(theta0)
        tilt=try(computeTiltDirect(b=b.3/N,g=gr,THRESH=0.00001,theta0=theta0, use.binomial=FALSE), silent=TRUE)
        if(!(class(tilt)=="try-error")) break
        numtries = numtries+1
    }    
    if(numtries==MAX.TRIES) stop("Newton Raphson could not determine correct tilt.")
    th = tilt$theta
    psi = tilt$psi
    psidot=tilt$psidot
    psidotdot = tilt$psidotdot
    
    mu = computeLocalIncrementMean(g=g,gdot=gdot,theta=tilt$theta,p0=p0)
    
    integrand<-function(t){
#        nuFunction(sqrt(2*N*mu)/sqrt(m*t*(1-t)))^2/(t^2*(1-t))
        (1-t)*nuFunction(sqrt(2*N*mu)/sqrt(m*t))^2/(t^2)
    }
    
    integral=integrate(integrand,lower=m0/m,upper=m1/m)
    
    cat("th=",th,", mu=",mu,", integral=", integral$value,"\n")
    
    # break the full expression below into two parts, to compare with pvalueSumChisqNew in the case p0=1:
    # exp(-N*(th*psidot-psi))*(2*pi*N*psidotdot)^(-1/2)*(th^(-1))*N^2*mu^2*integral$value
    part1=exp(-N*(th*psidot-psi))*(2*pi*N*psidotdot)^(-1/2)
    part2=(th^(-1))*N^2*mu^2*integral$value
    
    pval = part1*part2
    pval
}



pvalueSumChisqNew<-function(b,m,m0,m1,N){
    
    CbN.2 = 1-N/b
    
    integrand<-function(u){
        nuFunction(sqrt(b)*CbN.2/(sqrt(T)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral.2=integrate(integrand,lower=m0/m,upper=m1/m)
    
    part1=2*dchisq(b,N)
    part2=.25*b^2*CbN.2^3*integral.2$value
    
    pval=part1*part2
    
    pval
    
    #.5*b^2*CbN.2^3*dchisq(b,N)*integral.2$value
    
}

pvalueSumChisq<-function(b.normed,T,delta,T0,N){

    b2=sqrt(b.normed*sqrt(2*N)+N)
    
    delta1 = min(1, T0/T)
    CbN.1 = 1-(N-1)/b2^2
    
    integrand<-function(u){
        vu(b2*CbN.1/(sqrt(T)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral.1 =  integrate(integrand,lower=delta,upper=delta1)
    
    pval=2^(-1)*b2^4*(b2^(-1)*2^(-1))*dchi(b2,N)*CbN.1^3*integral.1$value
    pval
}

dchi<-function(y,N){
    fy <-(1-N/2)*log(2) + (N-1)*log(y) - y^2/2 - lgamma(N/2)
    exp(fy)
}

pmarg.sumweightedchisq<-function(b,ALPHA,N){
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0);
    gnormed<-function(u){(g(u)-g.moments$psidot)/sqrt(g.moments$psidotdot)}
    theta0=sqrt(g.moments$psidotdot)/2   # need theta/sqrt(psidotdot) < 1/2 for psi(theta)<infty.
    b0=b/sqrt(N)
    THRESH=0.0001
    marg = computeTiltDirect(b0,gnormed,THRESH,theta0)
    thetabN = marg$theta*sqrt(N)
    pmarg = exp(-thetabN*b + N*marg$psi)/sqrt(2*pi*marg$psidotdot)
    pmarg
}

computeMomentsBinomial<-function(g,theta,m){
# m is the window size, assumes that u~binomial(m,1/2).

    u = (seq(0,m,1)-m/2)/sqrt(m*.5^2)
    p = (0.5^m)*choose(m,seq(0,m,1))
  
    psidottop =  sum( g(u)*exp(theta*g(u))*p)
    psidotbot = sum(exp(theta*g(u))*p)
    psidot = psidottop/psidotbot
    EgUsq = sum(g(u)^2*exp(theta*g(u))*p)
    psidotdot = (EgUsq*psidotbot - psidottop^2)/(psidotbot^2)
    psi = log(psidotbot)

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}


computeMoments<-function(g,theta){


    INTLIM.THRESH=0.001
    psidottop.int <-function(u){ g(u)*exp(theta*g(u)-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=psidottop.int(INTLIM)
        if (!is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if (temp<INTLIM.THRESH) break
    }
    psidottop =  integrate(psidottop.int,lower=-INTLIM,upper=INTLIM)
    psidottop = psidottop$value/sqrt(2*pi)

    psidotbot.int =  function(u){ exp(theta*g(u)-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
    psidot = psidottop/psidotbot

    EgUsq.int<-function(u){ g(u)^2*exp(theta*g(u))*exp(-(u)^2/2)/psidotbot }    
    
    # us = seq(-INTLIM,INTLIM,.1)
    # plot(us,EgUsq.int(us))
    
    for( INTLIM in 10:50 ){
        temp=EgUsq.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH) break
    }
    
    EgUsq =  integrate(EgUsq.int,lower=-INTLIM,upper=INTLIM)
    EgUsq = EgUsq$value/sqrt(2*pi)
    
    if(EgUsq<0) stop("Error, overshot theta: E[g^2 e^(th*g-psi(th))]-E[g e^(th*g)]^2 < 0")

    psidotdot = EgUsq - (psidottop/psidotbot)^2;
    psi = log(psidotbot)

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}


computeTiltDirect <-function(b,gr,THRESH,theta0, use.binomial=FALSE,binomial.winsize=NA) {
#    Computes, via Newton Raphson, the value of theta such that
#    E_{theta}[g(U)]= b, for U~N(0,1),
#    and given function handle g.
#    E_{theta} is defined as the expectation under the directly
#    tilted measure for g(U), and not for U.
#    
#    Used in importance sampling of:
#    P(max_t sum_1^N g(U_ti) / sqrt(N) > b'), where b'/sqrt(N) = b.

#    cat("Trying Newton Raphson with theta0=",theta0,"\n")
    
    theta = theta0 # if start theta at 0, then dh(theta)=0 at the first step for g(u)=u^2.
    prevtheta=Inf
    prevprevtheta=Inf
    thetarec = theta
    
    while( abs(theta-prevtheta)>THRESH && abs(theta-prevprevtheta)>THRESH ){
        
        if(use.binomial && binomial.winsize>0) {
            g.moments<-computeMomentsBinomial(gr,theta,binomial.winsize)
        } else { 
            g.moments<-computeMoments(gr,theta)
        }
        htheta=g.moments$psidot-b
        dhtheta = g.moments$psidotdot

        # cat("theta=", theta," htheta=",htheta,".\n",sep="")
        thetarec= c(thetarec, theta)
        prevprevtheta=prevtheta
        prevtheta=theta
        theta = prevtheta - htheta/dhtheta
    }

    theta=prevtheta
    psi = g.moments$psi

    list(theta=theta,psi=psi,psidot=g.moments$psidot,psidotdot=g.moments$psidotdot)
}


vu<-function(x,maxn=1000,do.approx=(abs(x)<0.2)){
#    Evaluates the function vu(x) defined in (4.35) and (4.37) in 
#    Siegmund (1985) Sequential Analysis. 
#    
#    If x<0 then vu(x)= vu(-x).  x must be less than 0.
#    If approx = 1, uses the approximation:
#    v(x) = exp(-Px), where P=0.583.
#    maxn is the upper cap in evaluating the sum in the exponential.

    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

#    if(sum(is.na(do.approx)) > 0 && msscan.debug.trace){
#        #check for errors here
#        print(do.approx)
#        scan()
#        print(x)
#        scan()
#    }


    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583);
    }else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583);
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n=matrix(c(1:maxn),nrow=1,ncol=maxn)
        summands = pnorm(-0.5*matrix(x[notdo.approx.ix],nrow=length(notdo.approx.ix),ncol=1)%*%sqrt(n))/matrix(rep(n,length(notdo.approx.ix)),ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum");
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm);
    }

    vux
}

getCutoffEpidemicChangeLocationKnown<-function(pval){
    qchisq(1-pval,df=1)
}


computeBeta <- function(ALPHA){
    w<-function(u){
         exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    integrand<-function(u){
        u^2*w(u)^2*(2*u^4*w(u)*(1-w(u))+5*u^2*w(u)-3*u^2-2)*dchi(u,1)
    }
    
#    us=seq(0,10,0.1)
#    ys=integrand(us)
#    plot(us,ys)
    
    numerator=integrate(integrand,lower=0,upper=10)
    
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0)
    numerator$value/(2*g.moments$psidotdot)
}



intmode <- function(x){
# returns the mode of x, assuming it consists of only
# positive integers.
  xt<-tabulate(x)
  xmode<-which(xt == max(xt))
  xmode
}



msscan.old<-function(y,win,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,SAMPLE.PVAL.CUTOFF=0.0001,
                    SAMPLE.LOW.THRESH = 0.1, SAMPLE.HIGH.THRESH=0.5,
                    SAMPLE.CUTOFF.METHOD=4, verbose=TRUE,
                    WCHISQ.CUTOFF=NA, SAMPLE.THRESHOLD=NA,
                    MIN.SNPS = 2, Z=NULL){

                                        # --------------------------------------------------------------------
                                        # msscan
                                        #  
                                        # This function evolves from msseg4.
                                        #
                                        # --------------------------------------------------------------------
                                        # DESCRIPTION:
                                        # 
                                        # Multi-sample segmentation using the sum-of-chi-square statistic.
                                        # Does not recursively recompute Chisquare after each cut because
                                        # assumes that the aberrations are small and sparse and thus the global mean
                                        # is representative of the baseline mean.
                                        # 
                                        # INPUTS:
                                        # -------
                                        # y: N x T matrix of N samples at T locations.
                                        # win: maximum allowed window size.
                                        # ALPHA: weight.
                                        # GLOBAL_PVAL_CUTOFF: cutoff of p-value for global segmentation.
                                        # SAMPLE_CUTOFF_METHOD: sample cutoff is absolute intensity (1) or pvalue (2) or both (3,4)?
                                        # MIN.SNPS: Minimum length (in # snps) of CNV.
                                        #
                                        # OUTPUTS:
                                        # --------
                                        # yhat: the smoothed version of y
                                        # chpts: a M x 2 matrix the segments from yhat that are nonzero, i.e. aberrations.
                                        # chisq: a M x N matrix, where M is number of segments, with
                                        #          the chisq contribution of each sample for that segment.
                                        # Z: the random field of sum-of-chisq statistic.
                                        # ----------------------------------------------------------------------
    
    N=dim(y)[1]
    T=dim(y)[2]
    DELTA = 1/T
    
    if(is.na(WCHISQ.CUTOFF)){
        WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win,N,ALPHA)
    }
    
    if(SAMPLE.CUTOFF.METHOD==2 || SAMPLE.CUTOFF.METHOD==3 || SAMPLE.CUTOFF.METHOD==4){
        if(is.na(SAMPLE.THRESHOLD)){
            SAMPLE.THRESHOLD = getCutoffEpidemicChangeLocationKnown(SAMPLE.PVAL.CUTOFF)
        }
    }
    
    if(verbose){
        cat("Running msscan: ",N," samples, ",T," SNPs, ALPHA=",ALPHA,"\n", sep="")
        cat("Global p-value cutoff: ",GLOBAL.PVAL.CUTOFF,", corresponding Chisq cutoff: ",WCHISQ.CUTOFF, ".\n", sep="");
    }
    
#    if(msscan.debug.trace){
#         print("about to call ComputeZ.c")
#         print(paste("dim(y) = ", nrow(y), ncol(y)) )
#         print(paste("ALPHA = ", ALPHA) )
#         scan()
#    }

    if(is.null(Z)){
        Z = ComputeZ.C(y, T, win, ALPHA)
        Z[,1:MIN.SNPS-1] = 0
#        if(msscan.debug.trace){
#            print("about to call computeMoments")
#            scan()
#        }
        # use this gnorm to compute moments, which assumes u is a gaussian.
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    }
       
    Z.passinds = which(Z>WCHISQ.CUTOFF,arr.ind=TRUE)

    if(length(Z.passinds) == 0){
         y.hat <- y
         y.hat[,] <- mean(y)
         return( list(yhat=yhat, chpts=matrix(nrow=0, ncol=2), sampleswithsegment=NULL,chisq=NULL,Z=NULL)  )
    }

    Z.pass = Z[Z.passinds]
    Z.pass.order =order(Z.pass,decreasing=TRUE)
    Z.pass = Z.pass[Z.pass.order]
    Z.passinds = Z.passinds[Z.pass.order,]
    
    chpts = cbind(Z.passinds[,1]+1,Z.passinds[,1]+Z.passinds[,2])
    chpts.keep = Remove.Overlap.C(chpts)
    chpts = chpts[chpts.keep,]    

    if(length(chpts)==2) chpts=matrix(chpts,nrow=1)
    
#    if(msscan.debug.trace){
#         print("about to call Chisq.Contrib.C")
#         scan()
#    }

    chisq <- Chisq.Contrib.C(y, T, chpts) # chisq[i,j] holds contribution of sample j to segment i.

    yhat = matrix(0,nrow=N,ncol=T)
    sampleswithsegment=matrix(0,nrow=N,ncol=nrow(chpts))
    throw= rep(0,nrow(chpts))

    for(i in 1:nrow(chpts)){
        
        if(SAMPLE.CUTOFF.METHOD==1){
            sampleswithsegment[,i] = abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH    
        }
        if(SAMPLE.CUTOFF.METHOD==2){
            sampleswithsegment[,i] = chisq[i,]>SAMPLE.THRESHOLD    
        }
        if(SAMPLE.CUTOFF.METHOD==3){
            sampleswithsegment[,i] = (abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH) | (chisq[i,]>SAMPLE.THRESHOLD)
        }
        if(SAMPLE.CUTOFF.METHOD==4){
            sampleswithsegment[,i] = ((abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.LOW.THRESH) & (chisq[i,]>SAMPLE.THRESHOLD)) | (abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH)
        }


        m2 <- y[which(sampleswithsegment[,i]==1),chpts[i,1]:chpts[i,2]]
        if( sum(sampleswithsegment[,i]==1) == 1 || length(chpts[i,1]:chpts[i,2]) == 1){
            m2 <- matrix(m2, nrow=sum(sampleswithsegment[,i]==1), ncol=length(chpts[i,1]:chpts[i,2])  )
        }
        yhat[which(sampleswithsegment[,i]==1),chpts[i,1]:chpts[i,2]] = apply(m2, 1, mean)
        if( sum(sampleswithsegment[,i])==0) throw[i] = 1
        
    }
    
    # Throw away those chpts i where sum(sampleswithsegment[,i])=0.
    chpts = chpts[which(!throw),]
    sampleswithsegment = sampleswithsegment[,which(!throw)]
    chisq = chisq[which(!throw),]
    
    list(yhat=yhat, chpts=chpts, sampleswithsegment=sampleswithsegment,chisq=chisq,Z=Z)
}


msscan.mbic<-function(y,y.var,win,ALPHA=0,verbose=TRUE,MIN.SNPS = 2,MAX.CHPTS=100,MAX.SKIPPED=MAX.CHPTS,do.plots=TRUE, PRIOR.CARRIER.PROB=0.1,PRIOR.SEGCOUNT=1,COMMON.PHAT=FALSE){
# y is TxN, T is number of SNPs, N is number of samples.  
    
    y = t(y)
        
    N=dim(y)[1]
    T=dim(y)[2]
    if(verbose){
        cat("Running msscan.mbic: ",N," samples, ",T," SNPs, ALPHA=",ALPHA,"\n", sep="")
    }
    

    # ---- Compute Z matrix.  Z[t,k]: changed segment starting at t+1, ending at t+k. #
    Z = ComputeZ.C(y, T, win, ALPHA)
    if(MIN.SNPS>1) Z[,1:MIN.SNPS-1] = 0
    g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
    g2.moments=computeMoments(g2,0);
    Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
#    if(do.plots){    
#        par(mfrow=c(2,1))
#        heatmap(t(y))
#        image.plot(Z)
#    }
    Z.vec = matrix(data=t(Z),nrow=1)
    Z.vec.ind = matrix(data=c(1:prod(dim(Z))), nrow=nrow(Z),ncol=ncol(Z),byrow=TRUE)
            # check:  these two should be the same.
            #  Z.vec[Z.vec.ind[5,20:30]]
            #  Z[5,20:30]
    Z.vec.ord = order(Z.vec, decreasing=TRUE)
            # check: 
            # plot(Z.vec[Z.vec.ord])
    segbag = matrix(nrow=0,ncol=2)
    carriers = matrix(nrow=N,ncol=0)
    bic = rep(0,MAX.CHPTS)
    bic.all = matrix(nrow=N,ncol=MAX.CHPTS,data=0)
    sum.deltasq = rep(0,MAX.CHPTS)
    sum.deltasqwin = rep(0,MAX.CHPTS)
    Z.ind = 1
    Z.len = length(Z.vec)
    nskipped=0

    while(TRUE){
        # (1) Find the change-point with next largest Z that's not deleted in Z.vec.
        while(Z.ind<Z.len){
            if(Z.vec[Z.vec.ord[Z.ind]] != 0) break
            Z.ind = Z.ind+1
        }
        if(Z.ind >= Z.len || Z.vec[Z.vec.ord[Z.ind]]<0) break
        st = ceiling(Z.vec.ord[Z.ind]/win)
        w = Z.vec.ord[Z.ind]%%win
        if(w == 0) w = win
        ed = st + w
        st = st+1  # changed segment starts at st and ends at ed.
        if(verbose){
            cat("Considering segment: (",st,",",ed,")... ")
#            heatmap(t(y),loc=c((st-10):(ed+10)))
#            image.plot(Z[(st-10):(ed+10),])
        }
        
        chisq <- Chisq.Contrib.R(y, T, matrix(c(st-1,ed),nrow=1,ncol=2)) # chisq[i,j] holds contribution of sample j to segment i.
        samp.ord = order(chisq,decreasing=TRUE)  # the order that samples enter into bag.
        
        if(nrow(segbag)==0){
            bic.update = update.bic.baseline(t(y),segbag,carriers,rep(0,0),rep(0,0),c(st,ed),samp.ord,win,y.var, PRIOR.SEGCOUNT=5, do.plots=do.plots)
        } else {
            bic.update = update.bic.baseline(t(y),segbag,carriers,sum.deltasq[1:nrow(segbag)],
                                sum.deltasqwin[1:nrow(segbag)],c(st,ed),samp.ord,win,y.var,PRIOR.CARRIER.PROB=0.1,PRIOR.SEGCOUNT=1, do.plots=do.plots)
        }
        if(COMMON.PHAT){
            mbic.loc = bic.update$mbic.one.phat
        } else {    
            mbic.loc = bic.update$mbic.sep.phat
        }
        sum.deltasq.loc = bic.update$sum.deltasq.loc
        sum.deltasqwin.loc = bic.update$sum.deltasqwin.loc
        nhat = which.max(mbic.loc)
        
        if(mbic.loc[nhat]<0){
            # this change-point should not be added.
            # thus, overlaps with it should not be removed.
            if(verbose){
                cat("max bic over all subsets is negative.\n")
            }
            nskipped = nskipped+1;
        } else {
            # add change-point, update bic.
            bic.all[,nrow(segbag)+1] = mbic.loc
            bic[nrow(segbag)+1] = mbic.loc[nhat]
            sum.deltasq[nrow(segbag)+1] = sum.deltasq.loc[nhat]
            sum.deltasqwin[nrow(segbag)+1] = sum.deltasqwin.loc[nhat]
            
            # then, update segbag, carriers.
            segbag = rbind(segbag, c(st,ed))
            this.carriers= rep(0,N)
            this.carriers[samp.ord[1:nhat]] = 1
            carriers = cbind(carriers,this.carriers)
            
            # finally, remove overlaps from consideration.
            for(i in c(max((st-win+1),2):ed)){
                overlapping.w = max(1,st-i+1):win
                Z.vec[Z.vec.ind[i-1,overlapping.w]] = 0  # Because Z records segments in terms of  [st+1, win] 
            }
            if(verbose){
                cat("bic maxed at ",nhat," carriers, value = ",mbic.loc[nhat],".\n")
            }
        }
        if(verbose){
                cat("   ",nrow(segbag)," change-points in bag, ",nskipped," skipped.\n")
        }
        if(nrow(segbag)>=MAX.CHPTS  || nskipped>=MAX.SKIPPED) break
        Z.ind = Z.ind+1
    }
   
    # Choose the final number of change-points by maximizing MBIC, and populate yhat.
    
    yhat = matrix(nrow=nrow(y),ncol=ncol(y),data=apply(y,1,median)) 
    mhat = which.max(bic)
    if(bic[mhat]>0){
        # at least some CNVs are found, populate yhat.
        final.cnvs = segbag[1:mhat,]
        final.carriers = carriers[,1:mhat]
        for(i in 1:mhat){
            w=final.cnvs[i,2]-final.cnvs[i,1]+1
            if(w>1){
                if(sum(final.carriers[,i])>1){
                    yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=apply(y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2],drop=FALSE],1,mean))
                } else {
                    yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=mean(y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]))
                }
            } else {
                yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]])          
            }
        }
    }
    
    if(do.plots){
        par(mfrow=c(2,1))
        plot(bic, type="b")
        lines(rep(mhat,2) ,c(min(bic),max(bic)), col="black")
        heatmap(t(carriers), main=paste("Carriers of CNVs in bag, first",mhat," selected."))
    }
    list(yhat=t(yhat), segbag=segbag, segbag.carriers=carriers, bic.all = bic.all, bic, bic=bic, chpts=final.cnvs, carriers=final.carriers)
}

update.bic.baseline<-function(this.y,segbag,carriers,this.sum.deltasq,this.sum.deltasqwin,newchpt,samp.ord,win,y.var, 
                                PRIOR.CARRIER.PROB=0.1, PRIOR.SEGCOUNT = 10, do.plots=TRUE){
#    # for debugging:
#    this.y=t(y)          ##### remember to back transpose later.
#    newchpt=c(st,ed)
#    this.sum.deltasqwin=sum.deltasqwin[1:nrow(segbag)]
#    this.sum.deltasq = sum.deltasq[1:nrow(segbag)]
    
    
    T = nrow(this.y)
    N = ncol(this.y)
    
    m = nrow(segbag)+1    
    w = newchpt[2]-newchpt[1]+1
    
    # Each of these are length N vectors.
    if(newchpt[2]-newchpt[1]+1 >= 2){
        delta.hat = (apply(this.y[newchpt[1]:newchpt[2],],2,sum) - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    } else {
     delta.hat = (this.y[newchpt[1]:newchpt[2],] - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    }
    
    effect.hat = delta.hat^2/y.var
    sum.deltasq.loc = cumsum(effect.hat[samp.ord])  
    sum.deltasqwin.loc = cumsum(effect.hat[samp.ord]*w*(1-w/T))
    
    PRIOR.Nm = PRIOR.SEGCOUNT*N
    PRIOR.Carriers = PRIOR.CARRIER.PROB*PRIOR.Nm
    
    # Computing the prior probability of J, assuming a single phat=P(carrier) over all segments.
    sumJi = sum(carriers) + c(1:N) 
    Jprior.one.phat = sumJi*log((sumJi+PRIOR.Carriers)/(N*m+PRIOR.Nm)) + (N*m-sumJi)*log((N*m+PRIOR.Nm-sumJi-PRIOR.Carriers)/(N*m+PRIOR.Nm))
    
    # Computing the prior probability of J, assuming a different phat=P(carrier) over all segments.
    J = apply(carriers,2,sum)
    Jprior.sep.phat= sum( J*log((J+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-J)*log((N+PRIOR.Nm-J-PRIOR.Carriers)/(N+PRIOR.Nm)))
    Jprior.sep.phat = Jprior.sep.phat + c(1:N)*log((c(1:N)+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-c(1:N))*log((N+PRIOR.Nm-c(1:N)-PRIOR.Carriers)/(N+PRIOR.Nm))
         
    tauprior = - (lchoose(T,m) + log(m))    

    maxlik = 0.5*(sum(this.sum.deltasqwin) + sum.deltasqwin.loc)
    
    # do not take log if sum.deltasq, sum.deltasqwin is ZERO.
    
    betaprior.terms =  -sumJi/2 -(sumJi/2)*log((sum(this.sum.deltasqwin)+sum.deltasqwin.loc)/sumJi) 
    
    randomwalk.terms = 2*m*(1.6-1.5) - sum(log(this.sum.deltasq)) - log(sum.deltasq.loc)
    

    mbic.one.phat = maxlik + tauprior + Jprior.one.phat + betaprior.terms + randomwalk.terms 
    mbic.sep.phat = maxlik + tauprior + Jprior.sep.phat + betaprior.terms + randomwalk.terms 

    
    par(mfrow=c(3,1))
    plot(sqrt(effect.hat[samp.ord]*w*(1-w/T)),ylab="Chi")
    plot(maxlik,type="b", ylim=c(min(mbic.sep.phat, na.rm=TRUE),max(maxlik)))
    points(mbic.one.phat,type="b",col="red")
    segments(which.max(mbic.one.phat),0,which.max(mbic.one.phat), max(maxlik),col="red",lwd=1)
    points(mbic.sep.phat,type="b",col="blue")
    segments(which.max(mbic.sep.phat),0,which.max(mbic.sep.phat), max(maxlik),col="blue",lwd=1)
    plot(Jprior.one.phat,type="b",col="red", ylim=c(min(c(Jprior.sep.phat,betaprior.terms),na.rm=TRUE),0), ylab="Penalties",
        main=paste("PRIOR.SEGCOUNT=",PRIOR.SEGCOUNT,", PRIOR.CARRIER.PROB=",PRIOR.CARRIER.PROB))
    points(Jprior.sep.phat,type="b",col="blue")
    points(betaprior.terms,type="b",col="darkgreen")
    points(randomwalk.terms,type="b",col="orange")
     
    list(mbic.sep.phat = mbic.sep.phat, mbic.one.phat = mbic.one.phat, sum.deltasq.loc = sum.deltasq.loc, sum.deltasqwin.loc = sum.deltasqwin.loc)
}




mscbs.old <- function(y,win,MIN.SNPs=3,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,use.bic=TRUE, MAX.CHPTS=NA,WCHISQ.CUTOFF=NA,plots=TRUE){

                                        # Arguments:
                                        # ----------
                                        # y:  must be TxN, with each COLUMN being a sample.
                                        # win: The maximum size of each split, for higher speed and coarser results, reduce win.
                                        # MIN.SNPs: The minimum number of SNPs that must reside in each segmented region.
                                        # ALPHA: weight.
                                        # GLOBAL.PVAL.CUTOFF: Cutoff for when to stop segmenting.
                                        # MAX.CHPTS: Maximum number of change-points desired.
                                        # WCHISQ.CUTOFF: A pre-computed cutoff.
                                        # plots: Whether to plot the fitted segmentation after each recursion -- useful for diagnostics.


                                        # For debugging:
                                        #  win=300; Z=NULL; MIN.SNPs=3; ALPHA=0; MAX.CHPTS=NA; WCHISQ.CUTOFF=NA; plots=TRUE; GLOBAL.PVAL.CUTOFF=0.0001

  N=dim(y)[2]
  T=dim(y)[1]
  DELTA = MIN.SNPs/T
  
  if(is.na(MAX.CHPTS)) MAX.CHPTS=floor(T/MIN.SNPs)
  
  if(is.na(WCHISQ.CUTOFF)){
    WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win,N,ALPHA)
    cat("MSCBS: weighted chisquare cutoff = ",WCHISQ.CUTOFF,"\n")
    
  }

  zlim=c(-3,3)

  S = apply(y,2,cumsum)
  y.var<-compute.var(y)
  yhat = matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
  y.r = y-yhat
  if(plots) image.plot(c(1:T),c(1:N),yhat,zlim=zlim,xlab="position",ylab="sample",main="Fitted");
  
  chpts = c(1,T)
  bestZ =compute.max.Z.C(y.r,win,y.var,ALPHA,MIN.SNPs)
  best.subchpt= matrix(bestZ$bestchpt,ncol=1,nrow=2)
  best.Z = bestZ$bestZ


  splitnum=0
  chpt.hist = vector("list",MAX.CHPTS)
  bic.terms = vector("list",MAX.CHPTS)
  bic = rep(0,MAX.CHPTS)
  yhat.rec = vector("list",MAX.CHPTS)
     
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

    # Classify samples and update yhat.
    y.r.classify = mscbs.classify(y.r[chpts[max.region]:chpts[(max.region+1)],] , newchpt-chpts[max.region]+1 , y.var , CHISQ.PVAL.THRESH=0.001,)
    y.r.hat = matrix(0,nrow=T, ncol=N)
    y.r.hat[chpts[max.region]:chpts[(max.region+1)],] = y.r.classify$yhat
    yhat = yhat + y.r.hat
    y.r = y.r - y.r.hat
    if(plots) image.plot(c(1:T),c(1:N),yhat,zlim=zlim, xlab="position",ylab="sample",main="Fitted");
    
    
    if(!is.na(newchpt[2])){  # The added change consistes of two change-points.
      cat("Split ",splitnum,": ",newchpt[1],", ",newchpt[2],", Z-score = ",max.Z,".\n",sep="")
      y.r.L = y.r[chpts[max.region]:(newchpt[1]-1),]
      y.r.M =  y.r[newchpt[1]:(newchpt[2]-1),]
      y.r.R = y.r[newchpt[2]:chpts[max.region+1],]
      
      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.M = compute.max.Z.C(y.r.M,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
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
      
      
      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
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
    chpt.hist[[splitnum]] = list(chpts=chpts,max.region=max.region,
               newchpt=newchpt,max.Z=max.Z,carriers=y.r.classify$carriers)
    best.Z = c(leftpart.Z, best.Z.new, rightpart.Z)    
    best.subchpt = cbind(leftpart,best.subchpt.new,rightpart)
    chpts = c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])
    
    yhat.rec[[splitnum]] = yhat
    bic.terms[[splitnum]] = mbic.ms(S=S,chpts=chpts,sigma=sqrt(y.var))
    bic[splitnum] = bic.terms[[splitnum]]$bic
    
    
                                        # # For debugging: 
                                        #    cat(chpts)
                                        #    cat(rbind(best.Z, best.subchpt))
    
  }

  chpt.hist = chpt.hist[1:splitnum]
  bic.terms = bic.terms[1:splitnum]
  bic = bic[1:splitnum]
  term1 = rep(0,splitnum)
  term2 = rep(0,splitnum)
  term3 = rep(0,splitnum)
  for(i in 1:splitnum){
    term1[i] = bic.terms[[i]]$term1
    term2[i] = bic.terms[[i]]$term2
    term3[i] = bic.terms[[i]]$term3
  }
  
  
  
  if(splitnum>0){
        kstar = which.max(bic)
        chpts.final = c(chpt.hist[[kstar]]$chpts[1:chpt.hist[[kstar]]$max.region],
                        chpt.hist[[kstar]]$newchpt, 
                        chpt.hist[[kstar]]$chpts[(chpt.hist[[kstar]]$max.region+1):length(chpt.hist[[kstar]]$chpts)])  
        yhat.final = yhat.rec[[kstar]]
  } else {
    chpts.final=c(1,T)
  }
    


  if(length(chpts.final)>2) chpts.final = chpts.final[2:length(chpts.final)-1]
  else chpts.final = rep(0,0)

  list(chpt.hist=chpt.hist, chpts=chpts.final,yhat=yhat.final,bic=bic,bic.terms=bic.terms)
}


mbic.ms<-function(S,chpts,sigma, type=1, J=null, phat=NULL, include.all.terms=TRUE){
# chpts is [1,tau(1),...,tau(m),n].
# J[i,j] =1 if sample j changes at change-point i.

    N = ncol(S)
    m = length(chpts)-2

    if((is.null(J) && type != 1) || (type==1)){
        J = matrix(nrow=m+2,ncol=N,data=rep(1,N*(m+2)))
    }

#    if(type==1){
#        
#        # Assumes that all samples change at every change-point, 
#        # and does not include terms that are constant with respect to length
#        # of the sequence.
#    
#        U = matrix(data=0,ncol=N,nrow=m)
#        for(j in 1:m){
#            U[j,] = (S[chpts[j+1],] - (chpts[j+1]/chpts[j+2])*S[chpts[j+2],])/sqrt(chpts[j+1]*(1-chpts[j+1]/chpts[j+2]))
#            U[j,] = U[j,]/sigma
#        }
#        term1 = sum(U^2)/2
#        tdiff = chpts[2:(m+2)] - chpts[1:(m+1)]
#        term2 = -0.5*N*sum(log(tdiff))+0.5*N*log(nrow(S))
#        term3 = -m*lchoose(nrow(S),m)
#        term4=0
#        term5=0
#        term6=0
#        term7=0
#        
#        bic = term1+term2+term3
#        
#    }  else {

        if(m>0){
 
#            print(J)
        
            n.t = apply(J,2,next.tau)
            sum.J = apply(J,1,sum)
            
            # compute deltas.
            muhat = matrix(data=0,ncol=N,nrow=m+1)
            muhat[1,] = S[chpts[2],]/chpts[2]
            
            for(j in 1:m){
                muhat[j+1,]=(S[chpts[j+2],]-S[chpts[j+1],])/(chpts[j+2]-chpts[j+1])
            }
            deltahat = muhat[2:(m+1),] - muhat[1:m,]
            deltahat = deltahat*J[2:(m+1),]
            deltahat = deltahat/sigma
            if(is.null(dim(deltahat))) deltahat = matrix(data=deltahat,nrow=1,ncol=length(deltahat))
            deltahat.sq.sum = apply(deltahat^2,1,sum)
            deltahat.sq.sum = max(deltahat.sq.sum,1)    # lower bound at 1, so that log(...) would at least be non-negative.
            if(is.null(phat)) phat = sum(J[2:(m+1),])/(m*N)
                
            U = matrix(data=0,ncol=N,nrow=m)
            for(j in 1:m){
                U[j,] = (S[chpts[j+1],] - (chpts[j+1]/chpts[n.t[j+1,]])*diag(S[chpts[n.t[j+1,]],c(1:N)]))/sqrt(chpts[j+1]*(1-chpts[j+1]/chpts[n.t[j+1,]]))
                U[j,] = U[j,]/sigma
                U[j,] = J[j+1,]*U[j,] 
            }
            
            term1 = 0.5*sum(U^2)
#            term2 = -0.5*sum(sum.J[2:(m+1)]*log(sum(U^2)/sum(sum.J[2:(m+1)])))
            term2 = -0.5*sum(sum.J[2:(m+1)]*log(pmax(apply(U^2,1,sum), 0.00001)/pmax(sum.J[2:(m+1)],rep(1,m))))

            term3 =  -m*lchoose(nrow(S),m)
            term4 = -0.5*sum(sum.J[2:(m+1)])
            term5 = -sum(log(deltahat.sq.sum))
            term6 = 2.3*m-3*m/2
            if(phat == 0 || phat ==1){
                term7 = 0
            } else {
                term7 = m*N*(phat*log(phat)+(1-phat)*log(1-phat))
            }
        
            if(include.all.terms){
                bic = term1 + term2 + term3 + term4 + term5 + term6 + term7
            } else {
                bic = term1 + term2 + term3
            }
            
        } else {
            bic=0; term1=0; term2=0; term3=0; term4=0; term5=0; term6=0; term7=0;
        }
        
#    }
    
    
    list(bic=bic,term1=term1,term2=term2,term3=term3,term4=term4,term5=term5,term6=term6,term7=term7)
}


next.tau <- function(this.J){
# For a binary vector J, n.t[i] is the next time, after i,
# that J=1.
# Assumes that J is a vector with at least two elements,
# with J[1] and J[length(J)] all being 1.
    sel = which(this.J==1)
    n.t = rep(0,length(this.J))
    for(i in 2:length(sel)){
        n.t[sel[i-1]:(sel[i]-1)] = sel[i]
    }
    n.t[length(this.J)] = length(this.J) # this value doesn't matter, it will never be used.
    n.t
}


compute.var<-function(y, use.mean=TRUE){

    if(is.null(dim(y))){
         y = as.matrix(y,nrow=length(y),ncol=1)
    }

    N=dim(y)[2]
    T=dim(y)[1]
    y.diff<-y[1:(T-1),] - y[2:T,]
    
    if(!use.mean){
        if(N>1){
            y.var <- apply(y.diff^2,2,median)
        } else {
            y.var = median(y.diff^2)
        }
    } else {
        if(N>1){
            y.var <- apply(y.diff^2,2,mean)/2
        } else {
            y.var = mean(y.diff^2)/2
        }
    }
    
    y.var
}




mscbs.classify<-function(this.y,seg,y.var, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=NA, MIN.SUFF.ABSDIFF=NA){
#    cat("entered mscbs.classify!\n")
    
    N=dim(this.y)[2]
    T=dim(this.y)[1]

    if(is.na(MIN.REQ.ABSDIFF)){
      MIN.REQ.ABSDIFF = 0.5*sqrt(y.var)
    }
    if(is.na(MIN.SUFF.ABSDIFF)){
      MIN.SUFF.ABSDIFF = 5*sqrt(y.var)
    }
    
    if(is.na(seg[1]) && is.na(seg[2])){
      # Invalid change-point.
      sampleswithsegment = rep(0,N)
      this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
    } else {
      if(is.na(seg[2])){
        # Single change-point.
        # cat("Squeak!  ncol = ", ncol(this.y), " nrow=", nrow(this.y),"\n")
        sample.pval = computeZ.onechange.sample(t(this.y),seg[1],y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[1:seg[1],sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[1:seg[1],sampleswithsegment],nrow=seg[1]),2,mean),seg[1]),nrow=seg[1],byrow=TRUE)
        this.yhat[(seg[1]+1):T,sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[(seg[1]+1):T,sampleswithsegment],nrow=T-seg[1]),2,mean),T-seg[1]),nrow=T-seg[1],byrow=TRUE)        

      } else {
        # Two change-points, square wave change, see who has change in middle.
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[seg[1]:(seg[2]-1),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[seg[1]:(seg[2]-1),sampleswithsegment],nrow=seg[2]-seg[1]),2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
        this.yhat[c(1:(seg[1]-1),(seg[2]):T),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[c(1:seg[1],(seg[2]):T),sampleswithsegment], nrow=T-(seg[2]-seg[1])),2,mean),T-(seg[2]-seg[1])),nrow=T-(seg[2]-seg[1]),byrow=TRUE)

      }
    }
    
#     cat("exiting mscbs.classify!\n")
   
    
    list(carriers=sampleswithsegment, yhat=this.yhat)
}



compute.max.Z.C<-function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)
        
        this.Z<-ComputeZ.C(t(this.y), this.T, win, ALPHA)    
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        this.Z = (this.Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    

        ##*##
        this.Z[,1:(MIN.SNPs-1)] = 0
        this.Z[1:(MIN.SNPs-1),]=0
        for(i in 1:this.T-MIN.SNPs-1){
            maxwin = this.T-MIN.SNPs-i
            if(maxwin<win) this.Z[i, maxwin:win]=0
        }
        this.Z[(this.T-MIN.SNPs):this.T,]=0
        ##*##
        

        maxind = matrix.max(this.Z)
        bestchpt= c(maxind[1]+1,maxind[1]+maxind[2]+1)
        bestZ = this.Z[maxind[1],maxind[2]] # keep bestZ the best of the two-change-point scan.  If pruning out the left or right change-point made a big difference in bestZ, it would not have been pruned out anyway.    
        
        # Test the left change-point individually.
        if(bestZ!=0){
          pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
          pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval

         # Pruning of left and right change-points (Olshen and Venkatraman suggestion) happens here, even though we haven't yet tested whether the double change-point has significant p-value.  The point is, in case the double change-point has significant p-value, then we would need to test for the significance of the left and right change-points anyway.
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            bestchpt = c(bestchpt[2], NA)
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             bestchpt = c(bestchpt[1], NA)
            }
          }
        } else {
          bestchpt = c(NA,NA)
        }
    }
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}


feature.matrix.scan <- function(yhat,segments){

# Takes a NxT matrix yhat (output from msscan.r and a Mx2 list of segments,
# and returns the mean values in these segments.

    pos = segments[,1]
    feature.matrix = yhat[,pos]
    feature.matrix
}

feature.matrix.cbs <- function(yhat,chpts){
# Takes a TxN matrix yhat, and a vector of chpts, and
# returns the mean values between change-points.

  T = nrow(yhat)
  N= ncol(yhat)
  if(chpts[1] != 1) chpts = c(1,chpts,T)

  feature.matrix = matrix(0,nrow=N,ncol=length(chpts)-1)
  regions = matrix(nrow=length(chpts)-1, ncol=2,data=0)
  for(i in 1:(length(chpts)-1)){
    feature.matrix[,i] = yhat[chpts[i],]
    regions[i,] = c(chpts[i], chpts[i+1])
  }
  
  list(feature.matrix=feature.matrix, feature.regions=regions)
}



cnvtable<-function(chrom,pos,featuremat,segments,sampleswithsegment){
    nsegs = nrow(segments)
    
    NCOLS = 11
    cnvlist=matrix(0,nrow=0,ncol=NCOLS)
   
    for(i in 1:nsegs){
        st = segments[i,1]; ed = segments[i,2]
        carriers = which(sampleswithsegment[,i]==1)
        n.with.cnv = length(carriers)
        carriers.gain = carriers[which(featuremat[carriers,i]>0)]
        carriers.loss = carriers[which(featuremat[carriers,i]<0)]
        if(length(carriers.gain)>0) mean.gain = mean(featuremat[carriers.gain,i])
        else mean.gain=0
        if(length(carriers.loss)>0) mean.loss = mean(featuremat[carriers.loss,i])
        else mean.loss=0
        cnvlist = rbind(cnvlist,c(st,ed,ed-st+1,chrom[st],pos[st],pos[ed],length(carriers),length(carriers.gain),mean.gain,length(carriers.loss),mean.loss))
    }
    cnvlist=as.data.frame(cnvlist)
    names(cnvlist)<-c("Start index","End index","# SNPs","Chromosome","Start pos","End pos","# Carriers","# Gain","Mean gain","# Loss","Mean loss")
  
    cnvlist
}





norm.snp<-function(subdata, K=3){
# ----------------------------------------------------------- #
# Pre-processing functions.
# ----------------------------------------------------------- #

# subdata needs to have T rows, N columns for T SNPs, N samples.
# requires library(fields) to do image.plot.
# K is the number of clusters to use in k-means.


    subdata=t(subdata)
    nsnps = ncol(subdata)
    
    
    # First, center subdata so that each SNP has median 0, and each sample has median 0.
    meds=apply(subdata,2,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds,byrow=TRUE)
    temp=apply(subdata,1,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds, byrow=FALSE)
    
    
    
    # Perform k-means clustering on the samples.
    clust = kmeans(subdata, centers=K)

#    clust = pam(subdata,k=K) # can also try kmedoids.
#    par(mfrow=c(1,1))
#    barplot(clust$cluster)
#    
#    T=2000
#    loc =(-T+1):0
#
#    loc=loc+T
#    heatmap(t(subdata[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)
    
    # Within each cluster, normalize each SNP by its median.
    subdata.norm = matrix(0,nrow=nrow(subdata),ncol=nsnps)
    meds=matrix(0,nrow=K,ncol=nsnps)
    for(i in 1:K){
        inds = which(clust$cluster==i)
        # Subtract out median for each snp:
        meds[i,] = apply(subdata[inds,],2,median)
        meds.mat= matrix(meds[i,],nrow = length(inds),ncol=nsnps,byrow=TRUE)
        subdata.norm[inds,] = subdata[inds,] -meds.mat
        
#        par(mfrow=c(3,1))
#        loc=1:200
#        image(t(subdata[inds,loc]),zlim=c(-3,3))
#        image(t(meds.mat[,loc]),zlim=c(-3,3))
#        image(t(subdata.norm[inds,loc]),zlim=c(-3,3))
    }
    
    # Divide each SNP by its IQR.
    iqr = apply(subdata.norm, 2, IQR)
    iqr = iqr/median(iqr)
    iqr.mat = matrix(iqr, nrow=nrow(subdata),ncol=ncol(subdata), byrow=TRUE)
    subdata.norm = subdata.norm / iqr.mat
    
#    loc=1:200
#    par(mfrow=c(2,1))
#    image.plot(loc, 1:nrow(subdata),t(subdata[,loc]),xlab="SNP",ylab="Sample")
#    plot(loc,apply(subdata[,loc],2,median),type="l")
#    
#    par(mfrow=c(2,1))
#    colors=c("black","red","blue")
#    plot(meds[1,loc],type="l",col=colors[1])
#    for(i in 2:K) lines(meds[i,loc],col=colors[i])
#    plot(clust$centers[1,loc],type="l",col=colors[1])
#    for(i in 2:K) lines(clust$centers[i,loc],col=colors[i])
# Conclusion: not too much difference between median and mean.
    
#    snpid=125
#    hist(subdata[,snpid],20)



#    # look at the effects of normalization:
#       
#    T=2000
#    loc =(-T+1):0
#
#    loc=loc+T
#    par(mfrow=c(2,1))
#    heatmap(t(subdata[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)
#    heatmap(t(subdata.norm[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)


    list(normed=t(subdata.norm), cluster = clust$cluster, snp.medians = meds)
}


matrix.max <- function(Z){
    # If Z has NA values treat as -Inf
    na.inds = which(is.na(Z),arr.ind=TRUE)
    Z[na.inds]=-Inf

    row.max.col=max.col(Z, ties.method=c("random", "first", "last"))
    max.row = which.max(Z[cbind(1:nrow(Z),row.max.col)])
    c(max.row,row.max.col[max.row])
}

compute.max.Z<-function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
# This function computes the maximum of the weighted-sum-of-chisquare scan statistic.
# It normalizes the statistic to have mean 0 and variance 1.
    K=ncol(this.S)
    
    this.Z = ComputeZ.fromS.R(this.S, this.SST, this.imap, win, ALPHA, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
    g2.moments=computeMoments(g2,0);
    this.Z = (this.Z - K*g2.moments$psidot)/sqrt(g2.moments$psidotdot*K)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
 
}


compute.max.ProjectedZ<-function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){

# This function computes the maximum of the projected-chisquare scan statistic.
# It normalizes the statistic to have mean 0 and variance 1.

    K=ncol(this.S) # number of platforms.
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    
    this.Z = ComputeProjectedZ.fromS.R(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    this.Z = (this.Z - 1)/sqrt(2)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}


fcompute.max.ProjectedZ<-function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs, f=NULL){
# This function computes the maximum of the projected-chisquare scan statistic.
# It normalizes the statistic to have mean 0 and variance 1.


    K=ncol(this.S) # number of platforms.
    this.T = nrow(this.S)
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    win = min(win, this.T-1)
    
    temp = fscan.max(this.S, this.SST, this.imap=this.imap,MIN.SNPs=MIN.SNPs,rratio=rratio,f=f, use.Project.statistic=TRUE, verbose=FALSE)
    bestchpt = temp$seg
    bestZ = temp$maxZ
    
    list(bestchpt=bestchpt,bestZ=bestZ)
}


ComputeZ.fromS.anyg.R<-function(this.S,this.SST,this.imap,win,g,MIN.SNPs){
# This does the same thing as ComputeZ.fromS.R, except the user passes
# in any transformation g, instead of the parameter ALPHA.
#
# NOTE!!! g(.) works off U, not U^2, so the argument to g is assumed to be normal,
#         not chisquared!!

    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;
    for(i in 1:N){  # loop over samples.
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                SSb = nsnps*(diff1/nsnps-this.S[T,i]/totalsnps)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
                
                # 2010/02/21: need to change this part so that both left and the right segments have MIN.SNPs, not just the sum!!!
                # set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))  # previous code.
                nsnpsL = this.imap[1:(T-k),i] -this.imap[1,i]  # number of snps to the left of the first change-point.
                nsnpsR = totalsnps - this.imap[(k+1):T,i]+this.imap[1,i]
                set.to.zero = which(nsnps<MIN.SNPs | nsnpsL<MIN.SNPs | nsnpsR<MIN.SNPs)
                
                U[set.to.zero,k] = 0
        }
        U = sqrt(U)
        Z <- Z+g(U);
    }
    return(Z)

}


ComputeZ.fromS.R<-function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;
    for(i in 1:N){  # loop over samples.
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                SSb = nsnps*(diff1/nsnps-this.S[T,i]/totalsnps)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
                
                # 2010/02/21: need to change this part so that both left and the right segments have MIN.SNPs, not just the sum!!!
                # set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))  # previous code.
                nsnpsL = this.imap[1:(T-k),i] -this.imap[1,i]  # number of snps to the left of the first change-point.
                nsnpsR = totalsnps - this.imap[(k+1):T,i]+this.imap[1,i]
                set.to.zero = which(nsnps<MIN.SNPs | nsnpsL<MIN.SNPs | nsnpsR<MIN.SNPs)
                
                U[set.to.zero,k] = 0
        }
        Z <- Z+g(U);
    }
    return(Z)
}


ComputeZ.fromS.R.partial<-function(this.S,this.SST,this.imap,start.inds, end.inds,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          SSb = nsnps*(diff1/nsnps - this.S[T,]/totalsnps)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = (SSb/dfnum)/(SSw/dfden)
          
          # 2010/02/21: need to change this part so that both left and the right segments have MIN.SNPs, not just the sum!!!
          # set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          
          nsnpsL = this.imap[st,] - this.imap[1,] # number of snps to the left of the first change-point.
          nsnpsR = this.imap[T,] - this.imap[ed,]
          set.to.zero = which(nsnps<MIN.SNPs | nsnpsL<MIN.SNPs | nsnpsR<MIN.SNPs)

          U[set.to.zero] = 0
          Z[i,j] = sum(g(U))
        }
      }
    }

    Z
}


ComputeBYZ.fromS.R<-function(this.S,y.var,this.imap,win,delta,MIN.SNPs){
# Statistic proposed by Benny Yakir.
# delta = [ delta1, delta2], where one must be nonnegative, the other nonpositive.

    T = nrow(this.S) # number of snps.
    N = ncol(this.S) # number of samples.
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    l.delta.1=matrix(nrow=T, ncol=win, data=0);
    l.delta.2=matrix(nrow=T, ncol=win, data=0);
    Z=l.delta.1;
    for(i in 1:N){ 
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                l.delta.1[1:(T-k),k]  = delta[1]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
                l.delta.2[1:(T-k),k]  = delta[2]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
                              
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                l.delta.1[set.to.zero,k] = 0
                l.delta.2[set.to.zero,k] = 0
        }
        Z <- Z+pmax(l.delta.1,0)+pmax(l.delta.2,0)
    }
    return(Z)
}


ComputeBYZ.fromS.R.partial<-function(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          
          l.delta.1 = delta[1]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
          l.delta.2 = delta[2]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          l.delta.1[set.to.zero] = 0
          l.delta.2[set.to.zero] = 0
          Z[i,j] = sum(pmax(l.delta.1,0)+pmax(l.delta.2,0))
        }
      }
    }

    Z
}


ComputeProjectedZ.fromS.R<-function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
# Computes the projected Z statistic:
# Z = [delta ' X]^2 / delta'delta, where for a given s,t X[k] is the t statistic for testing for change in series k.

    T = nrow(this.S) # Number of SNPs.
    N = ncol(this.S) # Number of samples/platforms

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    weightU = U;
    sum.weight.sq = 0
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                temp<-diff1/nsnps-this.S[T,i]/totalsnps
                
                SSb = nsnps*(temp)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = sqrt(SSb/(SSw/dfden));
                weightU[1:(T-k),k] = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#                weightU[1:(T-k),k] = sign(temp)

                # 2010/02/21: need to change this part so that both left and the right segments have MIN.SNPs, not just the sum!!!
                # set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                nsnpsL = this.imap[1:(T-k),i]-this.imap[1,i] # number of snps to the left of the first change-point.
                nsnpsR = totalsnps - this.imap[(k+1):T,i]+this.imap[1,i]
                set.to.zero = which(nsnps<MIN.SNPs | nsnpsL<MIN.SNPs | nsnpsR<MIN.SNPs)
                U[set.to.zero,k] = 0
        }
        Z <- Z+rratio[i]*weightU*U;
        sum.weight.sq = sum.weight.sq + rratio[i]^2*weightU^2
    }
    Z<- Z^2/sum.weight.sq
    Z[which(is.na(Z), arr.ind=TRUE)]=0
    return(Z)
}


ComputeProjectedZ.fromS.R.partial<-function(this.S,this.SST,this.imap,start.inds, end.inds,rratio,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;


    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          temp = diff1/nsnps - this.S[T,]/totalsnps
          SSb = nsnps*(temp)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = sqrt((SSb/dfnum)/(SSw/dfden))
          
          weightU = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#          weightU = sign(temp)
          
          # 2010/02/21: need to change this part so that both left and the right segments have MIN.SNPs, not just the sum!!!
          # set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          nsnpsL = this.imap[st,] - this.imap[1,]# number of snps to the left of the first change-point.
          nsnpsR = this.imap[T,] - this.imap[ed,]
          set.to.zero = which(nsnps<MIN.SNPs | nsnpsL<MIN.SNPs | nsnpsR<MIN.SNPs)
          
          U[set.to.zero] = 0
          Z[i,j] = sum(U*weightU*rratio)/sqrt(sum(rratio^2*weightU^2))
        }
      }
    }
    
    Z<- Z^2
    Z[which(is.na(Z), arr.ind=TRUE)]=0
 
    Z
}




Chisq.Contrib.fromS.R <- function(S,SST,imap,chpts){
    T = nrow(S)
    N = ncol(S)
    dfnum <- 1;
    n.segs <- nrow(chpts)
    
    ret.m <- matrix(nrow=n.segs, ncol=N, data=0)
 
    totalsnps = imap[T,]-imap[1,]
    
    for(i in 1:n.segs){
        st <- chpts[i,1]
        ed <- chpts[i,2]
        w=ed-st
        
        for(j in 1:N){
            nsnps = imap[ed,j]-imap[st,j]        
            mn1 <- S[T,j]/totalsnps[j]
            SSb <- (S[ed,j]-S[st,j] -nsnps*mn1)^2 / (nsnps*(1-nsnps/totalsnps[j]));
            SSw = SST[j]-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/(totalsnps[j]-2));
        }

    }
    
    return(ret.m)
}


merge.pos<-function(pos,anchors=NULL){
    M = length(pos)
    pos.len = rep(0,M)
    for(i in 1:M) pos.len[i] = length(pos[[i]])
    
    if (is.null(anchors)) anchors=1:M
    
    nanchors = length(anchors)
    not.anchors = setdiff(c(1:M),anchors)
    
    totalinds = 0
    for(i in 1:nanchors){
        totalinds =totalinds + length(pos[[anchors[i]]]) 
    }
    
    merged.pos = rep(0,totalinds)
    index.map = matrix(0,nrow=totalinds,ncol=M)
    
    counter = rep(1,M) # counter[i] keeps track of where current merged.inds is in each of pos[i].
    next.pos = rep(0,nanchors)
    for(i in 1:nanchors) next.pos[i] = pos[[anchors[i]]][counter[anchors[i]]]
    
    for(k in 1:totalinds){
        min.pos = which.min(next.pos)
        merged.pos[k] = next.pos[min.pos]
        next.pos[min.pos] = pos[[anchors[min.pos]]][counter[anchors[min.pos]]+1]
        counter[anchors[min.pos]] = counter[anchors[min.pos]]+1
        for(j in not.anchors){
            while(TRUE) {
                if(counter[j] > pos.len[j]) break
                if(pos[[j]][counter[j]]>merged.pos[k]) break
                counter[j] = counter[j] + 1
            }
        }
        
        index.map[k,] = counter
    }

    # merged.pos is union of pos[[anchors]].
    # imap[k,i] = min_j { pos[[i]][j] > merged.pos[k] }   (note it is STRICTLY larger than)  
    # Thus to get the index of the largest value in pos[[i]] 
    # that is smaller than or equal to merged.pos[k], use imap[i,k]-1.
    # By this notation, imap[totalinds,] = pos.len + 1.

    list(merged.pos = merged.pos, imap = index.map)

}




fscan<-function(y,b,f=NULL,rho=0.99,MIN.SNPS=1, ALPHA=0, delta=c(1,-1), 
                use.BY.statistic=FALSE, throw.small.Z=FALSE, MIN.ABS.MEDIAN=0.2, 
                plots.pdf="fscanplots.pdf", verbose=FALSE){

# ------------------------------------------------------ #
# Sped up versions of scan by sequentially
# taking a denser anchor set.
# ------------------------------------------------------ #

  T=nrow(y)   # Number of SNPs
  N=ncol(y)   # Number of samples

  cat("Starting fscan with ",T," SNPs, ",N, " samples, threshold of ",b,".\n", sep="")

  if(is.null(f)){
  # This part doesn't work very well yet, initially, set f manually, e.g. f=c(0.01,0.1,1) works well.
    R= ceiling(log(log(T)) - log(-log(rho)))
    f = T^(-0.5^c(1:(R-1)))
    truncated.len = length(f)
    for(i in 2:length(f)){
      if(floor(1/f[i]) == floor(1/f[i-1])){
        truncated.len = i-1
        break
      }
    }
    f = c(f[1:truncated.len],1)
  }

  f = sort(f)
  if(f[length(f)] != 1) f = c(f,1)
  
  R = length(f)
  L = ceiling(f[2:R]/f[1:(R-1)])
  L = c(ceiling(T*f[1]),L)

  S = apply(y,2,cumsum)
  SST = apply(y,2,var)*(T-1)
  y.var <- compute.var(y)
 
  chpts = matrix(ncol=2,nrow=0)
  chpts.Z = matrix(nrow=1,ncol=0)

  g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
  g2.moments=computeMoments(g2,0);

  pdf(plots.pdf)
  
  for(r in 1:R){
    stepsize = floor(1/f[r])
    t = seq(1,T,stepsize) # t is the filtered anchor set.
    if(t[length(t)]<T) t = c(t,T) # always include the last datapoint in the set.

    this.S = S[t,]
    this.imap = matrix(rep(t,N),ncol=N,byrow=FALSE)

    # produce a diagnostic plot.
    plot(S[,1],xlab="SNP #", ylab="Cumulative sum of sample 1", main=paste("Round ",r,": Take 1 out of every ", stepsize," points.",sep=""))
    lines(t,this.S[,1],col="red")
    points(t,this.S[,1],col="red",pch=17,cex=2)

    # Refine previously found change-points using the denser anchor set.
    if(nrow(chpts)>0){
      for(i in 1:nrow(chpts)){
        ind.L = chpts[i,1] %/% stepsize
        ind.R = chpts[i,2] %/% stepsize
        
        check.win = f[r]/f[r-1]
        start.inds = c((ind.L - check.win):(ind.L+check.win))  
        end.inds = c((ind.R-check.win):(ind.R+check.win))
        
        if(use.BY.statistic){
            Z.part=ComputeBYZ.fromS.R.partial(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPS)
            Z.part = Z.part/N
        } else {
            Z.part=ComputeZ.fromS.R.partial(this.S,SST,this.imap,start.inds, end.inds,ALPHA,MIN.SNPS)
            Z.part = (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
        }

        #  image.plot(start.inds,end.inds,Z.part, xlab="Start of change - 1", ylab="End of change")

        maxind = matrix.max(Z.part)
        improved.cp = c(t[start.inds[maxind[1]]], t[end.inds[maxind[2]]])
        improved.Z = Z.part[maxind[1],maxind[2]]
        if(verbose) cat("Changepoints (",chpts[i,1],", ",chpts[i,2],") refined to (", improved.cp[1],", ",improved.cp[2],").  Z-score improved from ", chpts.Z[i]," to ",improved.Z,"\n", sep="")
        chpts[i,] = improved.cp
        chpts.Z[i] = improved.Z
      }
    }


    # Do scan with min window size MIN.SNPS and max window size L_i,    
    if(use.BY.statistic){
        Z = ComputeBYZ.fromS.R(this.S, y.var, this.imap, L[r], delta, MIN.SNPS)
        Z = Z/N
    } else {
        Z = ComputeZ.fromS.R(this.S, SST, this.imap, L[r], ALPHA, MIN.SNPS)
        Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    }


    # produce a diagnostic plot.
    image.plot(t,c(1:L[r]),Z,xlab="Start of change - 1", ylab="Window size", main=paste("Sum of chisquare (Z-score) for round ",r,sep=""))
    
    Z.passinds = which(Z>b,arr.ind=TRUE)

    if(length(Z.passinds) > 0){

      # Some locations passed the threshold!  Gather together all that passed,
      # remove the (abundant) overlaps, and merge with changepoints found in previous rounds.
      
      Z.pass = Z[Z.passinds]
      Z.pass.order =order(Z.pass,decreasing=TRUE)
      Z.pass = Z.pass[Z.pass.order]
      Z.passinds = Z.passinds[Z.pass.order,]
    
      newchpts = cbind(Z.passinds[,1],Z.passinds[,1]+Z.passinds[,2])
      if(length(newchpts) == 2) newchpts =matrix(nrow=1,ncol=2,data=newchpts)
      newchpts.keep = Remove.Overlap.C(newchpts)
      newchpts = newchpts[newchpts.keep,]
      newchpts = matrix(t[newchpts],ncol=2,byrow=FALSE)
      newZ = Z.pass[newchpts.keep]

      # This is specific to the scanning algorithm and is not necessary if CBS were used:
      # If, say, in a region 1-1000, the subregion 500-1000 has an amplification.  Then
      # without knowledge of a "baseline", 1-500 would look like a deletion (with very high
      # Z-score.  It would get called at this stage, and all subsequent focal changes within
      # 1-500 would be missed.  To guard against this, conduct a simple filter here.  If
      # none of the samples have a |median| > MIN.ABS.MEDIAN, then throw this change-point away.
      # In fact, if a large portion of the chromosome were changed, then one should not use the
      # simple scanning algorithm but instead use CBS.
      throw=rep(0,nrow(newchpts))
      for(i in 1:nrow(newchpts)){
        y.med = apply(y[newchpts[i,1]:newchpts[i,2],,drop=FALSE], 2,median)
        if(sum(abs(y.med)>MIN.ABS.MEDIAN)==0) throw[i] = 1
      }
      newchpts=newchpts[!throw,]
      if(length(newchpts)==2)newchpts=matrix(nrow=1,ncol=2,data=newchpts)
      newZ=newZ[!throw]
      ####
      
      
      if(verbose) {
        cat("Change-points found that are new to scan ",r,":\n")
        print(cbind(newchpts,newZ), digits=1)
      }

      chpts = rbind(chpts,newchpts)
      chpts.Z = c(chpts.Z,newZ)
      chpts.len = chpts[,2]-chpts[,1]

      # Check again the updated list for overlaps with CNVs found in previous rounds,
      # and if there are any, discard the one with
      # the smaller Z score (or the smaller length).  This can be sped up!!
      if(throw.small.Z){ 
        chpts.Z.order = order(chpts.Z, decreasing=TRUE)
      } else {
        chpts.Z.order = order(chpts.len, decreasing=TRUE)
      }
      chpts.Z = chpts.Z[chpts.Z.order]
      chpts = chpts[chpts.Z.order,]
      if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)
      chpts.keep = Remove.Overlap.C(chpts)
      chpts = chpts[chpts.keep,]
      if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)      
      chpts.Z = chpts.Z[chpts.keep]
      
    }

    if(verbose){
        cat("After merging, change-points found so far:\n")
        if(length(chpts.Z)>0) print(cbind(chpts, chpts.Z), digits=1)
        cat("\n")
    }
  }

  dev.off()
  
  list(chpts = chpts, chpts.Z = chpts.Z)
}


scan.classify <- function(y,segs, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
    nsamples=dim(y)[2]
    nsnps=dim(y)[1]
    y.var <- compute.var(y)
    yhat = matrix(nrow=nrow(y), ncol=ncol(y),data=0)
    if(length(segs)>2){ throw = rep(FALSE,nrow(segs))
    }else  {throw = FALSE}
    
    
    if(length(segs)==2){
        y.class = scan.classify.seg(y,segs,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)
        throw = sum(y.class)==0
    } else {
        y.class =matrix(ncol=nrow(segs),nrow=nsamples,data=FALSE)
        for(i in 1:nrow(segs)){
            seg=segs[i,]
            y.class[,i] <- scan.classify.seg(y,seg,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)  #  This is the step that takes the most time.
            throw[i] = sum(y.class[,i])==0
            if(sum(y.class[,i])==1){
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] = rep(mean(y[(seg[1]+1):seg[2],y.class[,i]]),seg[2]-seg[1])
                } else {
                    yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            } else {
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] <- matrix(rep(apply(y[(seg[1]+1):seg[2],y.class[,i]],2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
                } else {
                  yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            }
        }   
    }

    segs.f = segs[!throw,]
    y.class = y.class[,!throw]
    list(filtered.segs = segs.f, yhat = yhat, y.class = y.class)
}

scan.classify.seg<-function(this.y,seg,y.var,CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
        
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        if(seg[2]-seg[1]>1){ 
            sample.med.diff=apply(as.matrix(this.y[(seg[1]+1):seg[2],]),2,median) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        } else { 
            sample.med.diff=as.matrix(this.y[(seg[1]+1):seg[2],]) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        }
       
        sample.abs.med.diff = abs(sample.med.diff)
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = sample.abs.med.diff > MIN.SUFF.ABSDIFF
        cut1 = sample.abs.med.diff  < MIN.REQ.ABSDIFF

        sampleswithsegment = (pass1 | pass2) & (!cut1)
        sampleswithsegment

}


mscbs.get.fitted<-function(y,tauhat,carriers){
   m = length(tauhat)-2
    N = dim(y)[2]
   
    yhat = matrix(data=0,nrow=nrow(y),ncol=ncol(y))
    
  for(n in 1:N){ # loop over samples.
            st = 1
            while(TRUE){
                temp=which(carriers[(st+1):(m+2),n]==1)
                
                ed=temp[1]+st
                loc.st = tauhat[st]
                loc.ed = tauhat[ed]-1
                yhat[loc.st:loc.ed,n] = mean(y[loc.st:loc.ed,n])
                st = ed
                if(st==m+2){
                    yhat[tauhat[ed],n] =  yhat[tauhat[ed]-1,n] 
                    break
                }
            }
        }
    yhat
}


mscbs.classify.bic<-function(y,S=NULL,tauhat,y.var){
    m = length(tauhat)-2
    N = dim(y)[2]
    y.sd = sqrt(y.var)
    
#    dev.set(dev.prev())    
#    for(i in 1:length(tauhat)){
#        segments(tauhat[i],1,tauhat[i],N)
#    }
#    dev.set(dev.prev())
#    
    if(m==0){
        yhat = matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
        carriers = rep(0,0)
        phat = 0
    } else {
        Z = matrix(data=0,nrow=m,ncol=N)
        L = matrix(data=0,nrow=m,ncol=N)
        if(is.null(S)){
            S = apply(y,2,cumsum)
        } 
        
        for(j in 1:m){
            Z[j,] = S[tauhat[j+1],]-S[tauhat[j],] - ((tauhat[j+1] - tauhat[j])/(tauhat[j+2] - tauhat[j]))*(S[tauhat[j+2],]-S[tauhat[j],])
            L[j,] = (tauhat[j+1] - tauhat[j])*(1-(tauhat[j+1] - tauhat[j])/(tauhat[j+2] - tauhat[j]))
            Z[j,] = Z[j,]/sqrt(L[j,])
            Z[j,] = Z[j,]/y.sd
        }
        bX = matrix(Z/sqrt(L), nrow=prod(dim(Z)),ncol=1)
        n = matrix(L, nrow=prod(dim(Z)),ncol=1)
        
        M = min(length(bX),100)     # Added these two lines 11/23.  Because if p is not specified, BIC.means will attempt
                                    # to allocate an length(bX) x  length(bX) matrix, and cause memory problems.
        p=seq(1/(2*M),1-1/(2*M),1/(2*M))
        pbic = BIC.means(bX=bX,n=n, p=p, plots=FALSE)  # Use the empirical BIC to find the best fit.
        Jhat = rep(0,length(bX))
        Jhat[pbic$bic.J]=1
        Jhat = matrix(Jhat,nrow=m,ncol=N)
        phat = pbic$bic.phat
    
        yhat = matrix(data=0,nrow=nrow(y),ncol=ncol(y))
        Jhat.full = rbind(rep(1,N),Jhat,rep(1,N))
        for(n in 1:N){
            st = 1
            while(TRUE){
                temp=which(Jhat.full[(st+1):(m+2),n]==1)
                ed=temp[1]+st
#                loc.st = tauhat[st]
#                loc.ed = tauhat[ed]-1
                 loc.st = tauhat[st]+1  # changed 12/09/2009 from the above, when testing on leukemia data.
                loc.ed = tauhat[ed]
                yhat[loc.st:loc.ed,n] = mean(y[loc.st:loc.ed,n])
                st = ed
                if(st==m+2){
                    yhat[tauhat[ed],n] =  yhat[tauhat[ed]-1,n] 
                    break
                }
            }
        }
    
    }
    list(yhat=yhat,carriers=Jhat,phat=phat)
}


fcompute.max.Z<-function(this.y,win=Inf,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
# This speeds up compute.max.Z.C using filtering.
# Currently not using the y.var argument, but later
# may improve it to use that.
    

    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)

        this.S = apply(this.y,2,cumsum)
        this.SST = apply(this.y,2,var)*(this.T-1)
        
        # Find [t1, t2] that maximizes Z using filtered scan.
        temp = fscan.max(this.S,this.SST,this.imap=NULL,MIN.SNPs=MIN.SNPs,ALPHA=ALPHA)
        bestchpt= temp$seg
        bestZ = temp$maxZ
        
        # Test the left change-point individually.
        if(bestZ!=0){
        
            if(bestchpt[1]==1){ 
                pval.L = 1
            } else {
                pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
            }
            if(bestchpt[2]==this.T){
                pval.R =1
            } else {
                pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval
            }
         # Pruning of left and right change-points (Olshen and 
         # Venkatraman suggestion) happens here, even though we 
         # haven't yet tested whether the double change-point has 
         # significant p-value.  The point is, in case the double 
         # A change-point has significant p-value, then we would 
         # need to test for the significance of the left and 
         # right change-points anyway. 
         
         # 2010/02/21: Added throwLeft throwRight, because sometimes when the left fails, the right may also fail.
         # this is what happened in the Levinson 302 sample chrom 22 data, and it kept segmenting forever.
         # The left change-point had a low p-value whereas the right change-point is 1 away from this.T.
         throwLeft=FALSE; throwRight=FALSE
         
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            #bestchpt = c(bestchpt[2], NA)
            throwLeft = TRUE
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             # bestchpt = c(bestchpt[1], NA)
             throwRight=TRUE
            }
          }
          
          if(throwLeft && !throwRight){
            bestchpt = c(bestchpt[2], NA)
          }
          if(!throwLeft && throwRight){
            bestchpt = c(bestchpt[1], NA)
          }
          if(throwLeft && throwRight){
            bestchpt = c(NA,NA)
            bestZ = 0
          }
          
          
        } else {
          bestchpt = c(NA,NA)
        }
    }
    
    bestchpt = bestchpt+1  # Added 12/09/2009:  this function is only used in fmscbs, where it seems to be always off by 1.
    list(bestchpt=bestchpt,bestZ=bestZ)
}




fscan.max<-function(this.S, this.SST,  y.var=NULL, use.BY.statistic=FALSE, use.Project.statistic=FALSE, this.imap=NULL, 
                    f=NULL,MIN.SNPs=2,ALPHA=0,delta=c(1,-1), rratio=NULL, doplots=FALSE, verbose=FALSE){
# ------------------------------------------------------ #
# This does a filtered scan like fscan, but instead of
# returning all change-points with Z-score >b, it only
# returns the one pair [t1, t2] that maximizes Z-score.
# This is used in each round of fmscbs, via the call to
# fcompute.max.Z.  It should be faster than
# fscan, because by taking maxZ it does not have to
# get rid of overlaps.
# ------------------------------------------------------ #

  T=nrow(this.S)   # Number of SNPs
  N=ncol(this.S)   # Number of samples
  if(T<=MIN.SNPs){
    maxZ=NA
    seg=c(NA,NA)
  } else {

      if(is.null(f)){
        f.power = seq(min(-floor(log10(T))+1,0),0,1)
        f=10^f.power
      }
      
      if(T<1000) f=1 # even if a value for f is passed in, do not allow filtering for short sequences.
      
      f = sort(f)
      if(f[length(f)] != 1) f = c(f,1)
      R = length(f)
      L = ceiling(f[2:R]/f[1:(R-1)])
      L = c(ceiling(T*f[1]),L)
      
      chpts = matrix(ncol=2,nrow=0)
      chpts.Z = matrix(nrow=1,ncol=0)
    
      if(is.null(this.imap)){
        this.imap = matrix(rep(seq(1,T,1),N),ncol=N,byrow=FALSE)
      }
    
      g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
      g2.moments=computeMoments(g2,0);
      if(is.null(rratio)) rratio = rep(1,N)
      
      for(r in 1:R){
        stepsize = floor(1/f[r])
        
        t = seq(1,T,stepsize) # t is the filtered anchor set.
        if(t[length(t)]<T) t = c(t,T) # always include the last datapoint in the set.
    
        f.S = this.S[t,,drop=FALSE]
        f.imap = this.imap[t,,drop=FALSE]
    
        # produce a diagnostic plot.
        if(doplots){
          plot(this.S[,1],xlab="SNP #", ylab="Cumulative sum of sample 1", main=paste("Round ",r,": Take 1 out of every ", stepsize," points.",sep=""))
          lines(t,f.S[,1],col="red")
          points(t,f.S[,1],col="red",pch=17,cex=2)
        }
        
        # Refine previously found change-points using the denser anchor set.
        if(nrow(chpts)>0){
          for(i in 1:nrow(chpts)){
            ind.L = chpts[i,1] %/% stepsize
            ind.R = chpts[i,2] %/% stepsize
            
            check.win = f[r]/f[r-1]
            start.inds = c((ind.L - check.win):(ind.L+check.win))  
            end.inds = c((ind.R-check.win):(ind.R+check.win))
            
                  
            if(use.BY.statistic){
                Z.part=ComputeBYZ.fromS.R.partial(f.S,y.var,f.imap,start.inds, end.inds,delta,MIN.SNPS)
                Z.part = Z.part/N
            } else {
                if(use.Project.statistic){
                    Z.part=ComputeProjectedZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,rratio,MIN.SNPs)
                    Z.part = (Z.part -1)/sqrt(2)          
                } else {
                    Z.part=ComputeZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,ALPHA,MIN.SNPs)
                    Z.part = (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
                }
            }
    
            #  image.plot(start.inds,end.inds,Z.part, xlab="Start of change - 1", ylab="End of change")
    
            maxind = matrix.max(Z.part)
            improved.cp = c(t[start.inds[maxind[1]]]+1, t[end.inds[maxind[2]]])
            improved.Z = Z.part[maxind[1],maxind[2]]
            if(verbose) cat("fscan.max: Changepoints (",chpts[i,1],", ",chpts[i,2],") refined to (", improved.cp[1],", ",improved.cp[2],").  Z-score improved from ", chpts.Z[i]," to ",improved.Z,"\n", sep="")
            chpts[i,] = improved.cp
            chpts.Z[i] = improved.Z
          }
        }
        
        
        # Do scan with min window size MIN.SNPS and max window size L_i,    
        if(use.BY.statistic){
            Z = ComputeBYZ.fromS.R(f.S, y.var, f.imap, L[r], delta, MIN.SNPS)
            Z = Z/N
        } else {
            if(use.Project.statistic){
                Z = ComputeProjectedZ.fromS.R(f.S, this.SST, f.imap, L[r], rratio, MIN.SNPs)
                Z = (Z - 1)/sqrt(2)    
            } else {
                # Do scan with min window size MIN.SNPs and max window size L_i,
                Z = ComputeZ.fromS.R(f.S, this.SST, f.imap, L[r], ALPHA, MIN.SNPs)
                Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)      
            }
        }
        # produce a diagnostic plot.
        if(doplots) image.plot(t,c(1:L[r]),Z, xlab="Start of change - 1", ylab="Window size", main=paste("Sum of chisquare (Z-score) for round ",r,sep=""))
    
        maxind = matrix.max(Z)
        #  newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]])
        newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]-1])  # changed 10/22, if maxind[1]+maxind[2] might become length(t)+1...
        newZ = Z[maxind[1],maxind[2]]
        if(verbose){
          cat("fscan.max: Change-point found in round ",r,":\n")
          print(c(newchpt,newZ), digits=1)
        }
        chpts = rbind(chpts,newchpt)
        if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)      
        chpts.Z = c(chpts.Z,newZ)
      }
    
      # Take the maximum over the rounds, return chpt and Z score.
      max.ind = which.max(chpts.Z)
      maxZ = chpts.Z[max.ind]
      seg = chpts[max.ind,]
  }
  list(maxZ = maxZ, seg = seg)
}


scoreFit<-function(segments,reps, trios, fraction.overlap=0.8){
# segments[[i]] is an ncalls[i]x2 matrix of [start, end] points.
    falsecalls=vector("list",length(segments));
    
    numsegparams=0;
    
    for( i in 1:length(segments)){
         segdim=ncol(segments[[i]])
         if( !is.null(segdim) ) numsegparams = segdim
         
    }
    if(numsegparams==0){
        cat("No calls were made.\n")
        score=0
    } else {
    
        # Replicates
    
        for (i in 1:nrow(reps)){
            seg1=segments[[reps[i,1]]]
            seg2=segments[[reps[i,2]]]
            
            if(length(seg1)==numsegparams) {
                nseg1 = 1
                seg1 = matrix(nrow=1, ncol=2, data=seg1)
            } else {
                nseg1=nrow(seg1)
            }
            if(length(seg2)==numsegparams) {
                nseg2 = 1
                seg2 = matrix(nrow=1, ncol=2, data=seg2)
            } else {
                nseg2=nrow(seg2)
            }
                
            if(nseg1 > 0 && nseg2>0){
                matchedpairs=matrix(nrow=nseg1,ncol=nseg2, data=0)
                for (j in 1:nseg1){
                    for (k in 1:nseg2){
                        overlap=intersect(c(seg1[j,1]:seg1[j,2]),c(seg2[k,1]:seg2[k,2]))
                        overlaplen=length(overlap)
                        if( (overlaplen/(seg1[j,2]-seg1[j,1]+1) >=fraction.overlap) && (overlaplen/(seg2[k,2]-seg2[k,1]+1) >=fraction.overlap)){
                            matchedpairs[j,k] = 1
                        }
                    }
                }
                # falsecalls{i} contains the segments for sample i that are not found
                # in its replicate. 
                # falsecalls{i} = ( segment parameters , P(FP+FN) )
                # In the case of replicates, if that segment is not found in the replicate
                # sample, then P(FP or FN) = 1, it is either a false positive in this
                # sample, or a false negative in the other sample.
            
                falsecalls[[reps[i,1]]]=seg1[which(apply(matchedpairs,1,sum)==0),]  
                if(length(falsecalls[[reps[i,1]]])==numsegparams) falsecalls[[reps[i,1]]] = matrix(nrow=1, ncol=2, data=falsecalls[[reps[i,1]]] )
                falsecalls[[reps[i,1]]]=cbind( falsecalls[[reps[i,1]]], rep(1,sum(apply(matchedpairs,1,sum)==0)))
                
                falsecalls[[reps[i,2]]]=seg2[which(apply(matchedpairs,2,sum)==0),] 
                if(length(falsecalls[[reps[i,2]]])==numsegparams) falsecalls[[reps[i,2]]] = matrix(nrow=1, ncol=2, data=falsecalls[[reps[i,2]]] )
            
                falsecalls[[reps[i,2]]]=cbind( falsecalls[[reps[i,2]]], rep(1,sum(apply(matchedpairs,2,sum)==0)))
            } 
            
            if(nseg1==0 && nseg2>0){
                falsecalls[[reps[i,2]]]=seg2
                falsecalls[[reps[i,2]]]=cbind( falsecalls[[reps[i,2]]], rep(1,nseg2))
            }
            if(nseg1>0 && nseg2==0){
                falsecalls[[reps[i,1]]]=seg1
                falsecalls[[reps[i,1]]]=cbind( falsecalls[[reps[i,1]]], rep(1,nseg1))
    
            }
            
        }
        
        # Trios
        
        for (i in 1:nrow(trios)){
            child=segments[[trios[i,1]]]
            p1=segments[[trios[i,2]]]
            p2=segments[[trios[i,3]]]
            
            if(length(child)==numsegparams) {
                nchild = 1
                child = matrix(nrow=1, ncol=2, data=child)
            } else {
                nchild=nrow(child)
            }
            if(length(p1)==numsegparams) {
                np1 = 1
                p1 = matrix(nrow=1, ncol=2, data=p1)
            } else {
                np1=nrow(p1)
            }
            if(length(p2)==numsegparams) {
                np2 = 1
                p2 = matrix(nrow=1, ncol=2, data=p2)
            } else {
                np2=nrow(p2)
            }

            
            # ------ calculate the child-p1, child-p2, p1-p2 matrices of matched aberations -----
            
            if(nchild>0){
                
                
                cp1matchedpairs=matrix(nrow=nchild,ncol=np1, data=0)
                if(np1>0){
                    for( j in 1:nchild ){
                        for( k in 1:np1 ){
                            overlap=intersect(c(child[j,1]:child[j,2]),c(p1[k,1]:p1[k,2]))
                            overlaplen=length(overlap)
                            if( (overlaplen/(child[j,2]-child[j,1]+1) >=fraction.overlap) && (overlaplen/(p1[k,2]-p1[k,1]+1) >=fraction.overlap) ){
                                cp1matchedpairs[j,k] = 1
                            }
                        }
                    }
                }
                
                
                cp2matchedpairs=matrix(nrow=nchild,ncol=np2, data=0)
                if(np2>0){
                    for( j in 1:nchild){
                        for( k in 1:np2){
                            overlap=intersect(c(child[j,1]:child[j,2]),c(p2[k,1]:p2[k,2]))
                            overlaplen=length(overlap)
                            if( (overlaplen/(child[j,2]-child[j,1]+1) >=fraction.overlap) && (overlaplen/(p2[k,2]-p2[k,1]+1) >=fraction.overlap)){
                                cp2matchedpairs[j,k] = 1
                            }
                        }
                    }
                }
                
                # ----- calculate the false calls in child, p1, and p2 -----
                # for each of c, p1, p2, compare with the other two members in the
                # family.
                
                childfc = rep(0,0)
                for (j in 1:nchild){
                    if ((np1==0 || sum(cp1matchedpairs[j,])==0) && (np2==0 || sum(cp2matchedpairs[j,])==0)){
                        # child T p1 F p2 F => P(FP+FN)=1
                        childfc=rbind(childfc, c(child[j,], 1))
                    }
                }
                falsecalls[[trios[i,1]]]=childfc
            
            }
            
        }
        
        score=0
        for( i in 1:length(falsecalls)){
            if (length(falsecalls[[i]])>0){
                temp=falsecalls[[i]]
                score=score+sum(temp[,numsegparams+1])
            }
        }

    }

    list(num.inconsistent=score, falsecalls=falsecalls)

}



