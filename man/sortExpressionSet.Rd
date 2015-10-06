\name{sortExpressionSet}
\alias{sortExpressionSet}
\title{
Sort the order of samples for an ExpressionSet object
}
\description{
Sort the order of samples for an ExpressionSet object.
}
\usage{
sortExpressionSet(
    es, 
    varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat = c("\%m/\%d/\%Y", NA, NA)
)
}
\arguments{
  \item{es}{
An ExpressionSet.
}
  \item{varSort}{
A vector of phenotype variable names to be used to sort the samples
of \code{es}.
}
  \item{timeFormat}{
A vector of time format for the possible time variables in \code{varSort}.
The length of \code{timeFormat} should be the same as that
of \code{varSort}. For non-time variable, the corresponding time format
should be set to be equal to \code{NA}. Please refer to function  
\code{\link{strptime}} of the \code{base} package.
}
}
\value{
An ExpressionSet object with samples sorted based on the variables
indicated in \code{varSort}.
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
  

  es.sim2 = sortExpressionSet(
    es = es.sim, 
    varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat = c("\%m/\%d/\%Y", NA, NA)
  )
  
}
\keyword{ methods }
