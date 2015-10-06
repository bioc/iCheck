\name{glmWrapper}
\alias{glmWrapper}
\title{
Perform glm test for all gene probes
}
\description{
Perform glm test for all gene probes.
}
\usage{
glmWrapper(es, 
           formula = FEV1 ~ xi + age + gender, 
           pos.var.interest = 1,
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
  \item{formula}{
An object of class \code{formula}. 
The left handside of \code{~} is the response variable. 
Gene probe must be represented by the variable \code{xi}.
For example, \code{xi~age+gender} (gene probe is the response variable);
Or \code{FEV1~xi+age+gender} (gene probe is the predictor).
}
  \item{pos.var.interest}{ integer. Indicates which covariate
in the right-hand-size of \code{~} of \code{formula} is of
the interest.
pos.var.interest \eqn{= 0} means the intercept is of the interest.
If the covariate of the interest is an factor or interaction term 
with more than 2 levels, the smallest p-value will represent
the pvalue for the covariate of the interest.  
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
This function applies R function \code{glm}
for each gene probe.
}
\value{
A list with the following elements:
\item{n.sig }{Number of significant tests after p-value adjustment.}

\item{frame }{A data frame containing test results sorted according
to the ascending order of unadjusted p-values for the covariate 
of the interest. The data frame contains
7 columns: \code{probeIDs}, \code{geneSymbols} (gene symbols of the genes where the probes come from),
\code{chr} (numbers of chromosomes where the probes locate),
\code{stats} (z-value), \code{pval} (p-values of the tests for the covariate
of the interest),
\code{p.adj} (adjusted p-values), \code{pos} (row numbers of the probes in 
the expression data matrix).}

\item{statMat}{A matrix containing test statistics for all covariates
and for all probes. Rows are probes and columns are covariates.
The rows are ordered according
to the ascending order of unadjusted p-values for the 
covariate of the interest.
}

\item{pvalMat}{A matrix containing pvalues for all covariates
and for all probes. Rows are probes and columns are covariates.
The rows are ordered according
to the ascending order of unadjusted p-values for the covariate
of the interest.
}

\item{pval.quantile}{
Quantiles (minimum, 25%, median, 75%, and maximum) of p-values
for each covariate including intercept provided in the
input argument \code{formula}.

}

\item{frame.unsorted }{A data frame containing test results. 
The data frame contains
7 columns: \code{probeIDs}, \code{geneSymbols} 
(gene symbols of the genes where the probes come from),
\code{chr} (numbers of chromosomes where the probes locate),
\code{stats} (z-value for the covariate of the interest), 
\code{pval} (p-values of the tests for the covariate of the interest),
\code{p.adj} (adjusted p-values), \code{pos} (row numbers of the probes in 
the expression data matrix).}

\item{statMat.unsorted}{A matrix containing test statistics for all covariates
and for all probes. Rows are probes and columns are covariates.
}

\item{pvalMat.unsorted}{A matrix containing pvalues for all covariates
and for all probes. Rows are probes and columns are covariates.
}

\item{memGenes}{A numeric vector indicating the cluster membership 
of probes (unsorted). 
\code{memGenes[i]=1} if the \eqn{i}-th probe is significant (adjusted pvalue \eqn{<} \code{alpha}) with 
positive z-value for the covariate of the interest; 
\code{memGenes[i]=2} if the \eqn{i}-th probe is nonsignificant ; 
\code{memGenes[i]=3} if the \eqn{i}-th probe is significant with 
negative z-value for the covariate of the interest; 
}
\item{memGenes2}{
A numeric vector indicating the cluster membership 
of probes (unsorted). 
\code{memGenes2[i]=1} if the \eqn{i}-th probe is significant (adjusted pvalue \eqn{<} \code{alpha}). 
\code{memGenes2[i]=0} if the \eqn{i}-th probe is nonsignificant. 
}
\item{mu1}{
Mean expression levels for arrays for probe cluster 1 
(average taking across all probes with \code{memGenes} value equal to 1.
}
\item{mu2}{
Mean expression levels for arrays for probe cluster 2 
(average taking across all probes with \code{memGenes} value equal to 2.
}
\item{mu3}{
Mean expression levels for arrays for probe cluster 3 
(average taking across all probes with \code{memGenes} value equal to 3.
}

\item{resMat}{A matrix with \eqn{2p} columns, where \eqn{p} 
is the number of covariates (including intercept; for a nominal variable
with 3 levels say, there were 2 dummy covariates). 
The first \eqn{p} columns are p-values. The remaining \eqn{p}
columns are test statistics.
}
}
\author{
Weiliang Qiu <stwxq@channing.harvard.edu>,
        Brandon Guo <brandowonder@gmail.com>,
        Christopher Anderson <christopheranderson84@gmail.com>,
        Barbara Klanderman <BKLANDERMAN@partners.org>,
        Vincent Carey <stvjc@channing.harvard.edu>,
        Benjamin Raby <rebar@channing.harvard.edu>
}
\note{
If the covariate of the interest is a factor or interaction term
with more than 2 levels, then the p-value of the likelihood ratio test
might be more appropriate than the smallest p-value for 
the covariate of the interest.
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

   res.glm = glmWrapper(
  es = es.sim, 
  formula = xi~as.factor(memSubj), 
  pos.var.interest = 1,
  family = gaussian, 
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
