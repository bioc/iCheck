\name{genExprSet}
\alias{genExprSet}
\title{
Generate an ExpressionSet object
}
\description{
Generate a simple ExpressionSet object.
}
\usage{
genExprSet(
  ex, 
  pDat, 
  fDat = NULL, 
  annotation = "lumiHumanAll.db")
}
\arguments{
  \item{ex}{
A matrix of expression levels. Rows are gene probes
and columns are arrays. 
}
  \item{pDat}{
A data frame describing arrays. Rows are arrays and columns are variables
describing arrays. The row names of \code{pDat} must be the
same as the column of \code{ex}.
}
  \item{fDat}{
A data frame describing gene probes. Rows are gene probes and columns
are variables describing gene probes. The rownames of \code{fDat}
must be the same as that of \code{ex}.
}
  \item{annotation}{
character string. Indicating the annotation library 
(e.g. \code{lumiHumanAll.db} for the gene probes.
}
}
\value{
an ExpressionSet object.
}
\author{
Weiliang Qiu <stwxq@channing.harvard.edu>,
        Brandon Guo <brandowonder@gmail.com>,
        Christopher Anderson <christopheranderson84@gmail.com>,
        Barbara Klanderman <BKLANDERMAN@partners.org>,
        Vincent Carey <stvjc@channing.harvard.edu>,
        Benjamin Raby <rebar@channing.harvard.edu>

}
\keyword{ methods }
