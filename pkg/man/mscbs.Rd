\name{mscbs}
\alias{mscbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
    Multi-sample circular binary segmentation with p-value or BIC based model selection.
}
\description{
    Segments multiple aligned sequences and obtain a sparse representation as described in Zhang and Siegmund (2011).
}
\usage{
mscbs(y, use.bic = TRUE, MIN.SNPs = 3, ALPHA = 0, GLOBAL.PVAL.CUTOFF = 1e-04, MAX.CHPTS = NA, plots = TRUE, CARRIER.RETEST.THRESHOLD = 0.01, MIN.REQ.ABSDIFF = NA, MIN.SUFF.ABSDIFF = NA, CHISQ.PVAL.THRESH = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    \code{n} by \code{m} data matrix, for \code{n} locations and \code{m} samples.
}
  \item{use.bic}{
    TRUE if BIC should be used as stopping criterion.
}
  \item{MIN.SNPs}{
    The minimum number of SNPs in a segment.
}
  \item{ALPHA}{
    The weight used for the mixture model.
}
  \item{GLOBAL.PVAL.CUTOFF}{
    If BIC is not used, the p-value cutoff to be used for stopping the segmentation.
}
  \item{MAX.CHPTS}{
    If BIC is used, the maximum number of change-points to consider.
}
  \item{plots}{
  Show plots of progress of segmentation?  (If data set is huge and graphics is slow then this can cause significantly delay!)
}
  \item{CARRIER.RETEST.THRESHOLD}{
%%     ~~Describe \code{CARRIER.RETEST.THRESHOLD} here~~
}
  \item{MIN.REQ.ABSDIFF}{
%%     ~~Describe \code{MIN.REQ.ABSDIFF} here~~
}
  \item{MIN.SUFF.ABSDIFF}{
%%     ~~Describe \code{MIN.SUFF.ABSDIFF} here~~
}
  \item{CHISQ.PVAL.THRESH}{
%%     ~~Describe \code{CHISQ.PVAL.THRESH} here~~
}
}
\details{
    
    data(leukemia)
    seg=mscbs(y=leukemia,use.bic=TRUE, MIN.SNPs=5, MAX.CHPTS=30,plots=TRUE)
    par(mfrow=c(2,1))
    heatmap(pmax(pmin(y,1),-1), main="Raw", zlim=c(-1,1))
    heatmap(seg$yhat,zlim=c(-1,1),main="MSCBS-BIC Fitted")
    

}


\value{
   \item{yhat}{the fitted signal matrix for y}
   \item{chpts}{the change-point locations}
   \item{chpt.hist}{segmentation history: the order in which the change-points were called}
   \item{carriers}{which samples actually changed at each change-point?}
   \item{bic}{vector of BIC values of each split}
   \item{pve}{percent variance explained by the fit}
}
\references{
    Zhang, N.R., Siegmund, D.O., Ji, H., and Li, J., 2010, Detecting simultaneous change-points in multiple sequences. Biometrika, 97, 631--645.\cr
    Siegmund, D.O., Yakir, B. and Zhang, N.R., 2011, Detecting simultaneous variant intervals in aligned sequences.  Annals of Applied Statistics, in press.\cr
    Zhang, N.R. and Siegmund, D., 2011, Model Selection for High Dimensional, Multi-sequence Change-point Problems.  Submitted.\cr
}
\author{
Nancy R. Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    data(leukemia) 
    seg=mscbs(y=leukemia,use.bic=TRUE, MIN.SNPs=5, MAX.CHPTS=30,plots=TRUE) 
    par(mfrow=c(2,1)) 
    heatmap(pmax(pmin(leukemia,1),-1), main="Raw", zlim=c(-1,1)) 
    heatmap(seg$yhat,zlim=c(-1,1),main="MSCBS-BIC Fitted")  


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
