\name{pca3DPlot}
\alias{pca3DPlot}
\title{
Scatter plot of 3 specified principal components
}
\description{
Scatter plot of 3 specified principal components.
}
\usage{
pca3DPlot(pcaObj, 
          plot.dim = c(1,2, 3),
          labelVariable = "subjID", 
           hybName = "Hybridization_Name",
          outFileName = "test_pca_raw.pdf", 
          title = "Scatter plot of pcas", 
          plotOutPutFlag = FALSE, 
          mar = c(5, 4, 4, 2) + 0.1, 
          lwd = 1.5, 
          equalRange = TRUE, 
          xlab = NULL, 
          ylab = NULL, 
          zlab = NULL, 
          xlim = NULL, 
          ylim = NULL, 
          zlim = NULL, 
          cex.legend = 1.5, 
          cex = 1.5, 
          cex.lab = 1.5, 
          cex.axis = 1.5, 
          legendPosition = "topright", 
          ...)
}
\arguments{
  \item{pcaObj}{
An object returned by the function \code{pca} of the R package
\code{pcaMethods}.
}
  \item{plot.dim}{A vector of 3 positive-integer-value integer specifying
which 3 pcas will be plot.
}
  \item{labelVariable}{
The name of a column of the phenotype data matrix. The elements
of the column will replace the column names of the expression data matrix.   
}
 \item{hybName}{
character string. indicating the phenotype variable \code{Hybridization_Name}.
}
  \item{outFileName}{
Name of the figure file to be created.
}
  \item{title}{
Title of the scatter plot.
}
  \item{plotOutPutFlag}{
logical. \code{plotOutPutFlag=TRUE} indicates the plots will be output
to pdf format files. Otherwise, the plots will not be output
to external files.
}
  \item{mar}{
A numerical vector of the form 'c(bottom, left, top, right)'
          which gives the number of lines of margin to be specified on
          the four sides of the plot.  The default is 'c(5, 4, 4, 2) +
          0.1'. see \code{\link{par}}.
}
  \item{lwd}{
The line width, a _positive_ number, defaulting to '1'.
see \code{\link{par}}.
}
  \item{equalRange}{
logical. Indicating if the x-axis and y-axis have the same range.
}
  \item{xlab}{
Label of x axis.
}
  \item{ylab}{
Label of y axis.
}
  \item{zlab}{
Label of z axis.
}
  \item{xlim}{
Range of x axis.
}
  \item{ylim}{
Range of y axis.
}
  \item{zlim}{
Range of z axis.
}
  \item{cex.legend}{
Font size for legend.
}
  \item{cex}{
 numerical value giving the amount by which plotting text
          and symbols should be magnified relative to the default.
see \code{\link{par}}.
}
  \item{cex.lab}{
The magnification to be used for x and y labels relative to the current
setting of cex.
}
  \item{cex.axis}{
The magnification to be used for axis annotation relative to the current
setting of cex.

see \code{\link{par}}.
}
  \item{legendPosition}{
Position of legend.
Possible values are \dQuote{bottomright}, \dQuote{bottom}, 
\dQuote{bottomleft}, \dQuote{left}, \dQuote{topleft},
\dQuote{top}, \dQuote{topright}, \dQuote{right} and \dQuote{center}.
}
  \item{\dots}{
Arguments to be passed to \code{\link{plot}}.
}
}
\value{
A matrix of PCA scores. Each column corresponds to a principal component.
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

    pca.obj = getPCAFunc(es = es.sim, 
                     labelVariable = "subjID", 
                     hybName = "memSubj",
                     requireLog2 = FALSE,
                     corFlag = FALSE
)


pca3DPlot(pcaObj = pca.obj, 
          plot.dim = c(1,2,3),
          labelVariable = "subjID", 
          hybName = "memSubj",
          plotOutPutFlag = FALSE, 
          cex.legend = 0.5, 
          legendPosition = "topright") 
    

}
\keyword{ methods }
