\name{lkhrWrapper}
\alias{lkhrWrapper}
\title{
Perform glm test for all gene probes
}
\description{
Perform glm test for all gene probes.
}
\usage{
lkhrWrapper(es, 
           formulaReduced = FEV1 ~ xi + gender,
           formulaFull =    FEV1 ~ xi + age + gender,
           family = gaussian, 
           logit = FALSE, 
           pvalAdjMethod = "fdr", 
           alpha = 0.05, 
           probeID.var = "ProbeID", 
           gene.var = "Symbol", 
           chr.var = "Chromosome", 
           applier = lapply,
           verbose = TRUE) 
}
\arguments{
  \item{es}{
An LumiBatch object.
\code{fData(es)} should
contains information about probe ID, chromosome number and gene symbol.
}
  \item{formulaReduced}{
An object of class \code{formula}. 
Formula for reduced model.
The left handside of \code{~} is the response variable. 
Gene probe must be represented by the variable \code{xi}.
For example, \code{xi~gender} (gene probe is the response variable);
Or \code{FEV1~xi+gender} (gene probe is the predictor).
}
  \item{formulaFull}{
An object of class \code{formula}. 
Formula for Full model.
The left handside of \code{~} is the response variable. 
Gene probe must be represented by the variable \code{xi}.
For example, \code{xi~age+gender} (gene probe is the response variable);
Or \code{FEV1~xi+age+gender} (gene probe is the predictor).
}
  \item{family}{
By default is gaussian. refer to \code{\link{glm}}.
}
  \item{logit}{
logical. Indicate if the gene probes will be logit transformed.
For example, for DNA methylation data, one might want to logit transformation
for the beta-value (\eqn{methylated/(methylated+unmethylated)}).
}
  \item{pvalAdjMethod}{
One of p-value adjustment methods provided by
the R function \code{p.adjust} in R package \code{stats}:
\dQuote{holm}, \dQuote{hochberg}, \dQuote{hommel}, 
\dQuote{bonferroni}, \dQuote{BH}, \dQuote{BY}, \dQuote{fdr}, 
\dQuote{none}.
}
  \item{alpha}{
Significance level. A test is claimed to be significant
if the adjusted p-value \eqn{<} \code{alpha}.
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
  \item{applier}{
By default, it is lapply. If the library multicore is available,
can use mclapply to replace lappy.
}
  \item{verbose}{
logical. Determine if intermediate output need to be suppressed. By default
\code{verbose=TRUE},
intermediate output will be printed.
}

}
\details{
This function applies R functions \code{lrtest} in R package \code{lmtest}
and \code{glm} for each gene probe.
}
\value{
A list with the following elements:
\item{frame }{A data frame containing test results sorted according
to the ascending order of unadjusted p-values for the covariate 
of the interest. The data frame contains
8 columns: \code{probeIDs}, \code{geneSymbols} (gene symbols of the genes where the probes come from),
\code{chr} (numbers of chromosomes where the probes locate),
\code{Chisq} (chi square test statistic), 
\code{Df} (degree of freedom of the chisquare test statistic),
\code{pval} (p-values of the tests for the covariate
of the interest), 
\code{p.adj} (adjusted p-values), \code{pos} (row numbers of the probes in 
the expression data matrix).
The rows are ordered based on the descending order of chisquare test statistic.
}

\item{frame.unsorted }{A data frame containing test results. 
unordered \code{frame}.}
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

set.seed(1234567)
es.sim$age = rnorm(ncol(es.sim), mean=50, sd=5)
res.lkh = lkhrWrapper(
  es = es.sim, 
  formulaReduced = xi ~ memSubj,
  formulaFull =    xi ~ memSubj + age,
  family = gaussian(), 
  logit = FALSE, 
  pvalAdjMethod = "fdr", 
  alpha = 0.05, 
  probeID.var = "probe", 
  gene.var = "gene", 
  chr.var = "chr", 
  applier = lapply,
  verbose = TRUE) 

}
\keyword{ methods }
