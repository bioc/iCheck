\name{R2PlotFunc}
\alias{R2PlotFunc}
\title{
Draw heatmap of square of correlations among arrays
}
\description{
Draw heatmap of square of correlations among arrays.
}
\usage{
R2PlotFunc(
    es, 
    hybName = "Hybridization_Name",
    arrayType = c("all", "replicates", "GC"), 
    GCid = c("128115", "Hela", "Brain"),
    probs = seq(0, 1, 0.25), 
    col = gplots::greenred(75), 
    labelVariable = "subjID", 
    outFileName = "test_R2_raw.pdf", 
    title = "Raw Data R^2 Plot", 
    requireLog2 = FALSE, 
    plotOutPutFlag = FALSE, 
    las = 2, 
    keysize = 1, 
    margins = c(10, 10), 
    sortFlag = TRUE,
    varSort=c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat=c("\%m/\%d/\%Y", NA, NA),
    ...)
}
\arguments{
  \item{es}{
ExpressionSet object of QC probe profiles.
}
  \item{hybName}{
character string. indicating the phenotype variable \code{Hybridization_Name}.
}
  \item{arrayType}{
A character string indicating if the correlations
are calculated based on all arrays, arrays with replicates,
or genetic control arrays. 
}
  \item{GCid}{
A vector of character string. symbols for genetic control samples. The symbols
can be more than one.
}
  \item{probs}{
A vector of probabilities specify the quantiles of correlations
to be output.
}
  \item{col}{
colors used for the image. see the function \code{\link{heatmap.2}}
in R package \code{gplots}.
}
  \item{labelVariable}{
A character string indicating how to label the arrays.
}
  \item{outFileName}{
A character string.
The name of output file.
}
  \item{title}{
Title of the plot.
}
  \item{requireLog2}{
logical. \code{requiredLog2=TRUE} indicates probe expression levels
will be log2 transformed. Otherwise, no transformation will be performed.
}
  \item{plotOutPutFlag}{
logical. \code{plotOutPutFlag=TRUE} indicates the plots will be output
to pdf format files. Otherwise, the plots will not be output
to external files.
}
  \item{las}{
 'las' numeric in {0,1,2,3}; the style of axis labels.
          0 - always parallel to the axis,
          1 - always horizontal,
          2 - always perpendicular to the axis,
          or 3 - always vertical.

see \code{\link{par}}.
}
  \item{keysize}{
numeric value indicating the size of the key.
see the function \code{\link{heatmap.2}}
in R package \code{gplots}.
}
  \item{margins}{
numeric vector of length 2 containing the margins.
see the function \code{\link{heatmap.2}}
in R package \code{gplots}.
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
should be set to be equal to \code{NA}. The details of the time format
for time variable can be found in the R function \code{\link{strptime}}.
}
  \item{\dots}{
Arguments to be passed to \code{\link{heatmap.2}}.
}
}
\value{
A list with 3 elments. The first element \code{R2Mat} is 
the matrix of squared correlation.
The second element \code{R2vec} is the vector of 
the upper triangle of the matrix
of squared correlation (diagnoal elements are excluded).
The third element \code{R2vec.within.req} is the vector of 
within-replicate \eqn{R^2},
that is, any element in \code{R2vec.within.req} is the squared correlation
coefficient between two arrays/replicates for a subject.
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
  
    # draw heatmap for the first 5 subjects
    png(file="r2plot.png")
    R2PlotFunc(
      es = es.sim[, 1:5], 
      hybName = "memSubj",
      arrayType = c("all", "replicates", "GC"), 
      GCid = c("128115", "Hela", "Brain"),
      probs = seq(0, 1, 0.25), 
      col = gplots::greenred(75), 
      labelVariable = "subjID", 
      outFileName = "test_R2_raw.pdf", 
      title = "Raw Data R^2 Plot", 
      requireLog2 = FALSE, 
      plotOutPutFlag = FALSE, 
      las = 2, 
      keysize = 1, 
      margins = c(10, 10), 
      sortFlag = TRUE,
      varSort=c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
      timeFormat=c("\%m/\%d/\%Y", NA, NA))
    dev.off()
        
}
\keyword{ methods  }
