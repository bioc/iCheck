\name{densityPlots}
\alias{densityPlots}
\title{
Draw estimated density plots for all arrays
}
\description{
Draw estimated density plots for all arrays.
}
\usage{
densityPlots(
  es, 
  requireLog2 = TRUE, 
  myxlab = "expression level", 
  mymain = "density plots", 
  datExtrFunc = exprs, 
  fileFlag = FALSE, 
  fileFormat = "ps", 
  fileName = "densityPlots.ps")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
An \code{ExpressionSet} object that used to run the 
whole genome-wide tests.
}
  \item{requireLog2}{
logic. indicating if log2 transformation is required before estimating densities.
}
  \item{myxlab}{
character. indicating x-axis label.
}
  \item{mymain}{
character. indicating title of the plot.
}
  \item{datExtrFunc}{
name of the function  to extract genomic data. For
an \code{ExpressionSet} object, you should set
\code{datExtrFunc=exprs}; for a \code{MethyLumiSet} object,
you should set \code{datExtrFunc=betas}.
}
  \item{fileFlag}{
logic. indicating if plot should be saved to an external figure file.
}
  \item{fileFormat}{
character. indicating the figure file type. Possible values are \dQuote{ps}, \dQuote{pdf},
or \dQuote{jpeg}. All other values will produce \dQuote{png} file.
}
  \item{fileName}{
character. indicating figure file name (file extension should be specified). For example,
you set \code{fileFormat="pdf"}, then you can set \code{fileName="test.pdf"}, but not
\code{fileName="test"}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list object, the \eqn{i}-th element is the object returned by
function \code{density} for the \eqn{i}-th array.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Weiliang Qiu <stwxq@channing.harvard.edu>,
        Brandon Guo <brandowonder@gmail.com>,
        Christopher Anderson <christopheranderson84@gmail.com>,
        Barbara Klanderman <BKLANDERMAN@partners.org>,
        Vincent Carey <stvjc@channing.harvard.edu>,
        Benjamin Raby <rebar@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  # generate simulated data set from conditional normal distribution
  set.seed(1234567)
  es.sim = genSimData.BayesNormal(nCpGs = 100, 
    nCases = 20, nControls = 20,
    mu.n = -2, mu.c = 2,
    d0 = 20, s02 = 0.64, s02.c = 1.5, testPara = "var",
    outlierFlag = FALSE, 
    eps = 1.0e-3, applier = lapply) 
  print(es.sim)

  densityPlots(
    es = es.sim, 
    requireLog2 = FALSE, 
    myxlab = "expression level", 
    mymain = "density plots", 
    datExtrFunc = exprs, 
    fileFlag = FALSE, 
    fileFormat = "ps", 
    fileName = "densityPlots.ps")
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
