\name{LumiBatch2Table}
\alias{LumiBatch2Table}
\title{
Output slots (exprs, pData, fData) of an LumiBatch object into 3 text files
}
\description{
Output slots (exprs, pData, fData) of an LumiBatch object into 3 text files.
}
\usage{
LumiBatch2Table(
  es, 
  probeID.var="ProbeID",
  gene.var="Symbol",
  chr.var="Chromosome",
  sep = ",", 
  quote = FALSE,
  filePrefix = "test", 
  fileExt = "csv")
}
\arguments{
  \item{es}{
An LumiBatch object
}
  \item{probeID.var}{
character string. Name of the variable indicating probe ID in feature
data set.
}
  \item{gene.var}{
character string. Name of the variable indicating gene symbol in feature
data set.
}
  \item{chr.var}{
character string. Name of the variable indicating chromosome number in feature
data set.
}
  \item{sep}{
Field delimiter for the output text files
}
  \item{quote}{
  logical. Indicating if any character or factor. See also \code{\link{write.table}}.  
  }
  \item{filePrefix}{
Prefix of the names of the output files.
}
  \item{fileExt}{
File extension of the names of the output files.
}
}
\details{
Suppose \code{filePrefix="test"} and \code{fileExt=".csv"}.
Then, the file names of the 3 output files are:
\dQuote{test_exprs.csv}, \dQuote{test_pDat.csv}, 
and \dQuote{test_fDat.csv}, respectively.
}
\value{
None.
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

   LumiBatch2Table(
    es = es.sim, 
    probeID.var="probe",
    gene.var="gene",
    chr.var="chr",
    sep = ",", 
    quote = FALSE,
    filePrefix = "test", 
    fileExt = "csv")

}
\keyword{ methods }
