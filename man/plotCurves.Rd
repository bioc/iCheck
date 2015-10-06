\name{plotCurves}
\alias{plotCurves}
\title{
Plot trajectories of probe profiles across arrays
}
\description{
Plot trajectories of probe profiles across arrays
}
\usage{
plotCurves(
    dat, 
    curveNames, 
    fileName,
    plotOutPutFlag=FALSE,
    requireLog2 = FALSE, 
    cex = 1, 
    ylim = NULL, 
    xlab = "", 
    ylab = "intensity", 
    lwd = 3, 
    main = "Trajectory plot", 
    mar = c(10, 4, 4, 2) + 0.1,
    las = 2,
    cex.axis=1,
    ...)
}
\arguments{
  \item{dat}{
Numeric data matrix. Rows are probes; columns are arrays.
}
  \item{curveNames}{
Probe names.
}
  \item{fileName}{
file name of output figure.
}
  \item{plotOutPutFlag}{
logical. \code{plotOutPutFlag=TRUE} indicates the plots will be output
to pdf format files. Otherwise, the plots will not be output
to external files.
}

  \item{requireLog2}{
logical. \code{requiredLog2=TRUE} indicates probe expression levels
will be log2 transformed. Otherwise, no transformation will be performed.
}
  \item{cex}{
 numerical value giving the amount by which plotting text
          and symbols should be magnified relative to the default.
see \code{\link{par}}.
}
  \item{ylim}{
Range of y axis.
}
  \item{xlab}{
Label of x axis.
}
  \item{ylab}{
Label of y axis.
}
  \item{lwd}{
The line width, a _positive_ number, defaulting to '1'.
see \code{\link{par}}.
}
  \item{main}{
Main title of the plot.
}
  \item{mar}{
A numerical vector of the form 'c(bottom, left, top, right)'
          which gives the number of lines of margin to be specified on
          the four sides of the plot.  The default is 'c(5, 4, 4, 2) +
          0.1'. see \code{\link{par}}.
}
  \item{las}{
 'las' numeric in {0,1,2,3}; the style of axis labels.
          0 - always parallel to the axis,
          1 - always horizontal,
          2 - always perpendicular to the axis,
          or 3 - always vertical.

see \code{\link{par}}.
}
  \item{cex.axis}{
The magnification to be used for axis annotation relative to the current
setting of cex.

see \code{\link{par}}.

}
  \item{\dots}{
Arguments to be passed to \code{\link{plot}}.
}
}
\value{
no return value.
}
\author{
Weiliang Qiu <stwxq@channing.harvard.edu>,
        Brandon Guo <brandowonder@gmail.com>,
        Christopher Anderson <christopheranderson84@gmail.com>,
        Barbara Klanderman <BKLANDERMAN@partners.org>,
        Vincent Carey <stvjc@channing.harvard.edu>,
        Benjamin Raby <rebar@channing.harvard.edu>
}
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
    
  # plot trajectories of the first 5 genes
  plotCurves(
  dat = exprs(es.sim)[1:5,], 
  curveNames = featureNames(es.sim)[1:5], 
  plotOutPutFlag=FALSE,
  cex = 0.5,
  requireLog2 = FALSE) 

}
\keyword{ methods }
