\name{getPCAFunc}
\alias{getPCAFunc}
\title{
Get principal components of arrays
}
\description{
Get principal components of arrays.
}
\usage{
getPCAFunc(es, 
           labelVariable = "subjID", 
            hybName = "Hybridization_Name",
           requireLog2 = TRUE,
           corFlag = FALSE
)
}
\arguments{
  \item{es}{
An \code{ExpressionSet} object
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
  \item{corFlag}{
logical. Indicating if correlation matrix (\code{corFlag=TRUE}) or covariance (\code{corFlag=FALSE}) is used to obtain 
principal components.
}

}
\value{
A list with 3 elements:
  \item{es.s }{An \code{ExpressionSet} object with the arrays sorted according
to Batch_Run_Date, Chip_Barcode, and Chip_Address
 }
  \item{pcs }{An object returned by the function \code{prcomp} of the R
package \code{stats}. It contans the following components. \code{sdev} (the square roots of the eigenvalues 
of the covariance/correlation matrix); \code{rotation} (a matrix whose columns contain the eigenvectors);
  \code{x} (a matrix whose columns contain principal components); \code{center} (the centering used or \code{FALSE});
\code{scale} (the scale used or \code{FALSE})}
  \item{requireLog2 }{logical. The same value as the input \code{requireLog2}.}
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
}
\keyword{ methods }
