\name{msscan}
\alias{msscan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
   Multi-sample Scan for Change in Mean
}
\description{
    Scans aligned sequences from multiple samples for segments with increased or decreased mean value in a subset of the samples.
}
\usage{
msscan(y, pval=0.001, maxwin, ALPHA = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    \code{n} by \code{m} data matrix, for \code{n} locations and \code{m} samples.

}
  \item{pval}{
    The p-value threshold for the scan.
}
  \item{maxwin}{
    The maximum window size used in the scan.  If unspecified, the total length of the sequences will be used.
}
  \item{ALPHA}{
    The weight that should be used in transforming the chi-square statistic (see Siegmund, Yakir and Zhang (2011)).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{chpts }{chpts[i,1] and chpts[i,2] are the start and end of the i-th most significant segment.}
%%  \item{chpts.Z }{The "Z" statistic as defined in Siegmund, Yakir and Zhang (2011).  This is the square of the usual normal "z-score".}
%% ...
}
\references{
    Zhang, N.R., Siegmund, D.O., Ji, H., and Li, J., 2010, Detecting simultaneous change-points in multiple sequences. Biometrika, 97, 631--645. \cr
    Siegmund, D.O., Yakir, B. and Zhang, N.R., 2011, Detecting simultaneous variant intervals in aligned sequences.  Annals of Applied Statistics, in press. \cr
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
    data(sqcp)
    normed = norm.snp.pc(sqcp)  # this data has significant artifacts, which we can remove by pc normalization.
    sig = msscan(normed$normed,pval=0.0001,maxwin = nrow(normed$normed))
    par(mfrow=c(2,1))
    heatmap(normed$normed, zlim=c(-1.5,1.5), main="Normalized Data")
    sigline = rep(0,nrow(normed$normed))
    for(i in 1:length(sig$chpts.Z)){
        sigline[sig$chpts[i,1]:sig$chpts[i,2]] = sqrt(sig$chpts.Z[i])
    }
    plot(1:nrow(normed$normed),sigline, xlab="Probe order", ylab="z-score", type="l", xaxs="i", main="Z-scores of Detected Intervals")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Scan statistic}
\keyword{Change-point estimation}% __ONLY ONE__ keyword per line
