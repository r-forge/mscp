\name{norm.snp.pc}
\alias{norm.snp.pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Normalization for array based CGH data.
}
\description{
    Normalization by (1) removing low-rank latent factor by PCA, (2) optionally scale each SNP so that the inter-quartile ranges are the same (3) optionally centering each SNP and sample.
}
\usage{
norm.snp.pc(subdata, num.pc = 2, do.plots = TRUE, center.snp = TRUE, center.sample = TRUE, iqr.snp = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subdata}{
    \code{n} by \code{m} data matrix, for \code{n} locations and \code{m} samples.

}
  \item{num.pc}{
    number of latent factors to use.
}
  \item{do.plots}{
    plot the intermediate results.
}
  \item{center.snp}{
    should we center each snp after normalization?
}
  \item{center.sample}{
    should we center each sample after normalization?
}
  \item{iqr.snp}{
  should each SNP be scaled to have the same inter-quartile range?
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{normed }{The normalized matrix}
    \item{artifacts}{The artifact matrix}
    \item{pc}{The latent factors}
    \item{snp.medians}{The estimated medians of the SNPs}
}
\references{
    

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
