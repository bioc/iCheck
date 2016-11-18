\name{quantilePlot}
\alias{quantilePlot}
\title{
Plot trajectories of quantiles across arrays
}
\description{
Plot trajectories of quantiles across arrays.
}
\usage{
quantilePlot(
    dat, 
    fileName, 
    probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), 
    plotOutPutFlag = FALSE, 
    requireLog2 = FALSE, 
    sortFlag = TRUE,
    cex = 1, 
    ylim = NULL, 
    xlab = "", 
    ylab = "intensity", 
    lwd = 3, 
    main = "Trajectory plot of quantiles", 
    mar = c(15, 4, 4, 2) + 0.1, 
    las = 2, 
    cex.axis = 1)
}
\arguments{
  \item{dat}{
Expression data. Rows are gene probes; columns are arrays.
}
  \item{fileName}{
File name of output figure.
}
  \item{probs}{
quantiles (any real values between the interval \eqn{[0, 1]}).
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
  \item{sortFlag}{
logical. \code{sortFlag=TRUE} indicates arrays will be sorted by 
the ascending order of MAD (median absolute deviation)
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
Charater string. main title of the plot.
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
}
\value{
The quantile matrix with row quantiles and column array. 
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


   png(file="qplot.png")

     quantilePlot(
       dat = exprs(es.sim), 
       probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), 
       plotOutPutFlag = FALSE, 
       requireLog2 = FALSE, 
       sortFlag = TRUE)

   dev.off()
  
}
\keyword{ methods }
