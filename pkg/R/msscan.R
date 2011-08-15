msscan <-
function(y,pval=0.001,maxwin=NA,ALPHA=0,MIN.SNPS=1,...){
    if(is.na(maxwin)){ maxwin = nrow(y)}
    b = getCutoffMultisampleWeightedChisq(pval,nrow(y),1/nrow(y),maxwin,ncol(y),ALPHA)
    fscan(y,f=c(0.01,0.1,1), ALPHA = ALPHA,b=b, verbose=FALSE,MIN.SNPS=MIN.SNPS,...)
}

