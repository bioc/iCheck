\name{genSimData.BayesNormal}
\alias{genSimData.BayesNormal}
\title{
Generating simulated data set from conditional normal distributions
}
\description{
Generating simulated data set from conditional normal distributions.
}
\usage{
genSimData.BayesNormal(
  nCpGs, 
  nCases, 
  nControls,
  mu.n = -2,
  mu.c = 2,
  d0 = 20, 
  s02 = 0.64,
  s02.c = 1.5,
  testPara = "var", 
  outlierFlag = FALSE,
  eps = 0.001, 
  applier = lapply)
}
\arguments{
  \item{nCpGs}{
integer. Number of genes.
}
  \item{nCases}{
integer. Number of cases.
}
  \item{nControls}{
integer. Number of controls.
}
  \item{mu.n}{
numeric. mean of the conditional normal distribution for controls. See details.
}
  \item{mu.c}{
numeric. mean of the conditional normal distribution for cases. See details.
}
  \item{d0}{
integer. degree of freedom for scale-inverse chi squared distribution.
See details.
}
  \item{s02}{
numeric. scaling parameter for scale-inverse chi squared distribution
for controls. See details.
}
  \item{s02.c}{
numeric. scaling parameter for scale-inverse chi squared distribution
for cases. See details.
}
  \item{testPara}{
character string. indicating if the test is for testing equal mean,
equal variance, or both.
}
  \item{outlierFlag}{
logical. indicating if outliers would be generated. If \code{outlierFlag=TRUE},
then we followed Phipson and Oshlack's (2014) simulation studies to 
generate one outlier for each CpG site by replacing the DNA methylation 
level of one diseased subject by the maximum of the DNA methylation levels 
of all CpG sites.
}
  \item{eps}{
numeric. if \eqn{|mean0-mean1|<eps} then we regard \eqn{mean0=mean1}.
Similarly, 
if \eqn{|var0-var1|<eps} then we regard \eqn{var0=var1}.
\eqn{mean0} and \eqn{var0} are the mean and variance of
the chi squared distribution for controls.
\eqn{mean1} and \eqn{var1} are the mean and variance of
the chi squared distribution for cases.
}
  \item{applier}{
function name to do \code{apply} operation.
}
}
\details{
Based on Phipson and Oshlack's (2014) simulation algorithm.
 For each CpG site, variance of the DNA methylation was first sampled from an
   scaled inverse chi-squared distribution with degree of freedom 
\eqn{d_0} and scaling parameter \eqn{s_0^2}:
      \eqn{\sigma^2_i ~ scale-inv \chi^2(d_0, s_0^2)}.
  M value for each CpG was then sampled from a normal distribution
    with mean \eqn{\mu_n} and variance equal to the simulated variance 
\eqn{\sigma^2_i}.
For cases, the variance was first generated from
      \eqn{\sigma^2_{i,c} ~ scale-inv \chi^2(d_0, s_{0,c}^2)}.
  M value for each CpG was then sampled from a normal distribution
    with mean \eqn{\mu_c} and variance equal to the simulated variance 
\eqn{\sigma^2_{i,c}}.
}
\value{
An ExpressionSet object. The phenotype data of the ExpressionSet object
contains 2 columns: \code{arrayID} (array id) and {memSubj} (subject
membership, i.e., case (\code{memSubj=1}) or control (\code{memSubj=0})).
The feature data of the ExpressionSet object contains 4 elements:
\code{probe} (probe id), \code{gene} (psuedo gene symbol), \code{chr}
(psuedo chromosome number), and \code{memGenes}
(indicating if a gene is differentially expressed (when \code{testPara="mean"})
or indicating if a gene is differentially variable (when \code{testPara="var"}) ).

}
\references{
    Phipson B, Oshlack A. 
    DiffVar: A new method for detecting differential variability with application to methylation in cancer and aging. 
    \emph{Genome Biol} 2014; 15:465
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

}
\keyword{ method }
