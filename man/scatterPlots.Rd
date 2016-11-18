\name{scatterPlots}
\alias{scatterPlots}
\title{
Draw scatter plots for top results in whole-genome-wide analysis
}
\description{
Draw scatter plots for top results in whole-genome-wide analysis to
test for the association of probes to a continuous-type phenotype
variable.
}
\usage{
scatterPlots(
  resFrame, 
  es, 
  col.resFrame = c("probeIDs", "stats", "pval", "p.adj"), 
  var.pheno = "bmi", 
  outcomeFlag = FALSE,
  fitLineFlag = TRUE,
  var.probe = "TargetID", 
  var.gene = "Symbol", 
  var.chr = "Chr", 
  nTop = 20, 
  myylab = "expression level", 
  datExtrFunc = exprs, 
  fileFlag = FALSE, 
  fileFormat = "ps", 
  fileName = "scatterPlots.ps")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resFrame}{
A data frame stores testing results, which must contain columns
that indicate probe id, test statistic, p-value and optionally 
adjusted p-value.
}
  \item{es}{
An \code{ExpressionSet} object that used to run the 
whole genome-wide tests.
}
  \item{col.resFrame}{
A vector of characters indicating column names of \code{resFrame}
corresponding to probe id, test statistic, p-value and optionally
adjusted p-value.
}
  \item{var.pheno}{
character. the name of continuous-type phenotype variable that is used to
test the association of this variable to probes.
}
  \item{outcomeFlag}{
logic. indicating if \code{var.pheno} is the outcome variable in regression analysis.
}
  \item{fitLineFlag}{
logic. indicating if a fitted line \eqn{y=a+bx} should be plotted. If \code{outcomeFlag=TRUE},
then \eqn{y} is \code{var.pheno} and \eqn{x} is the top probe.
If \code{outcomeFlag=FALSE},
then \eqn{y} is the top probe and \eqn{x} is \code{var.pheno}. 
}
  \item{var.probe}{
character. the name of feature variable indicating probe id.
}
  \item{var.gene}{
character. the name of feature variable indicating gene symbol.
}
  \item{var.chr}{
character. the name of feature variable indicating chromosome number.
}
  \item{nTop}{
integer. indicating how many top tests will be used to draw
the scatter plot.
}
  \item{myylab}{
character. indicating y-axis label.
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
Value \code{0} will be returned if no error occurs.
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

  # generate phenotype age
  es.sim$age = rnorm(ncol(es.sim), mean=50, sd=5)

  res.limma = lmFitWrapper(
    es = es.sim, 
    formula = ~age, 
    pos.var.interest = 1,
    pvalAdjMethod = "fdr", 
    alpha = 0.05, 
    probeID.var = "probe", 
    gene.var = "gene", 
    chr.var = "chr", 
    verbose = TRUE)

  scatterPlots(
    resFrame=res.limma$frame, 
    es=es.sim, 
    col.resFrame = c("probeIDs", "stats", "pval"), 
    var.pheno = "age", 
    outcomeFlag = FALSE,
    fitLineFlag = TRUE,
    var.probe = "probe", 
    var.gene = "gene", 
    var.chr = "chr", 
    nTop = 20, 
    myylab = "expression level", 
    datExtrFunc = exprs, 
    fileFlag = FALSE, 
    fileFormat = "ps", 
    fileName = "scatterPlots.ps")
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
