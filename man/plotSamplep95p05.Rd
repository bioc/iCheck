\name{plotSamplep95p05}
\alias{plotSamplep95p05}
\title{
Plot trajectories of the ratio of 95th percentile
to 5th percentile of sample probe profiles across arrays
}
\description{
Plot trajectories of the ratio of 95th percentile
to 5th percentile of sample probe profiles across arrays.
}
\usage{
plotSamplep95p05(
    es, 
    labelVariable = "subjID", 
     hybName = "Hybridization_Name",
    requireLog2 = FALSE, 
    projectName = "test", 
    plotOutPutFlag = FALSE, 
    cex = 1, 
    ylim = NULL, 
    xlab = "", 
    ylab = "", 
    lwd = 1.5, 
    mar = c(10, 4, 4, 2) + 0.1, 
    las = 2, 
    cex.axis=1.5,
    title = "Trajectory of p95/p05",
    cex.legend = 1.5,
    cex.lab = 1.5,
    legendPosition = "topright",
    cut1 = 10,
    cut2 = 6,
    sortFlag = TRUE,
    varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat = c("\%m/\%d/\%Y", NA, NA),
    verbose = FALSE,
    ...)
}
\arguments{
  \item{es}{
ExpressionSet object of Sample probe profiles. 
}
  \item{labelVariable}{
A character string.
The name of a phenotype data variable use to
label the arrays in the boxplots. By default, 
\code{labelVariable = "subjID"} which is equivalent
to \code{labelVariable = "Hybridization_Name"}.
}
 \item{hybName}{
character string. indicating the phenotype variable \code{Hybridization_Name}.
}
  \item{requireLog2}{
logical. \code{requiredLog2=TRUE} indicates probe expression levels
will be log2 transformed. Otherwise, no transformation will be performed.
}
  \item{projectName}{
A character string. Name of the project. The plots will be saved  as pdf format files,
the names of which have the format
\code{projectName_probeName_traj_plot.pdf}.
}
  \item{plotOutPutFlag}{
logical. \code{plotOutPutFlag=TRUE} indicates the plots will be output
to pdf format files. Otherwise, the plots will not be output
to external files.
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
  \item{title}{
Figure title.
}
  \item{cex.legend}{
Font size of legend text.
}
  \item{cex.lab}{
The magnification to be used for x and y labels relative to the current
setting of cex.
}
  \item{legendPosition}{
Position of legend.
Possible values are \dQuote{bottomright}, \dQuote{bottom}, 
\dQuote{bottomleft}, \dQuote{left}, \dQuote{topleft},
\dQuote{top}, \dQuote{topright}, \dQuote{right} and \dQuote{center}.

}
  \item{cut1}{
  second horiztonal line setting the cutoff for the ratio \code{p95/p05}.
A ratio above this line indicates the corresponding array is good.
}
  \item{cut2}{
  second horiztonal line setting the cutoff for the ratio \code{p95/p05}.
A ratio below this line indicates the corresponding array is bad.

}
  \item{sortFlag}{
logical. Indicates if arrays need to be sorted according to 
\code{Batch_Run_Date}, \code{Chip_Barcode}, and \code{Chip_Address}.

}
  \item{varSort}{
A vector of phenotype variable names to be used to sort the samples
of \code{es}.
}
  \item{timeFormat}{
A vector of time format for the possible time variables in \code{varSort}.
The length of \code{timeFormat} should be the same as that
of \code{varSort}. For non-time variable, the corresponding time format
should be set to be equal to \code{NA}.
The details of the time format
for time variable can be found in the R function \code{\link{strptime}}.
}
  \item{verbose}{
logical. Determine if intermediate output need to be suppressed. By default
\code{verbose=FALSE},
intermediate output will not be printed.
}
  \item{\dots}{
Arguments to be passed to \code{\link{plot}}.
}
}
\details{
The trajectory of the ratio of 95% percentile
to 5% percentile of the expression levels of all QC probes will be ploted.
}
\value{
A list of 2 elements. The first element is the 2 x n matrix, where n
is the number of arrays. The first row of the matrix is the 5-th percentile
and the second row of the matrix is the 95-th percentile. 

The second element is the ratio of the 95-th percentile to the 5-th percentile.
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

  es.sim$Batch_Run_Date = 1:ncol(es.sim)
  es.sim$Chip_Barcode = 1:ncol(es.sim)
  es.sim$Chip_Address = 1:ncol(es.sim)
  

 plotSamplep95p05(
  es = es.sim, 
  labelVariable = "subjID", 
  hybName = "memSubj",
  requireLog2 = FALSE, 
  projectName = "test", 
  plotOutPutFlag = FALSE, 
  title = "Trajectory of p95/p05",
  cex.legend = 0.5,
  legendPosition = "topright",
  sortFlag = TRUE,
  varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
  timeFormat = c("\%m/\%d/\%Y", NA, NA),
  verbose = FALSE)

}
\keyword{ methods }
