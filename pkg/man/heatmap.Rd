\name{heatmap}
\alias{heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Heatmap Visualization for Array Copy Number Data
}
\description{
    Draws a heatmap of the data matrix with samples as rows and probes as columns.
}
\usage{
heatmap(y, chrom = NULL, loc = NULL, zlim = NULL, horizontal = TRUE, pos.string = "Probe order", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    n by m data matrix, with each column corresponding to a sample and each row corresponding to a probe location.
}
  \item{chrom}{
    the chromosome that each probe belongs to, a vector of length n.
}
  \item{loc}{
    The locations of the probes.
}
  \item{zlim}{
    The scale for the z (color) scale.
}
  \item{horizontal}{
    Should the heatmap be described horizontally (samples as rows) or vertically (samples as columns)?
}
  \item{pos.string}{
    What should the location axis be called?
}
  \item{\dots}{
    Other standard plotting parameters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
    heatmap(y=leukemia, zlim=c(-1.5,1.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
