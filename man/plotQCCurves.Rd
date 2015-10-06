\name{plotQCCurves}
\alias{plotQCCurves}
\title{
Plot trajectories of specific QC probes (e.g., biotin, cy3_hyb, housekeeping
gene probes, low stringency probes, etc.) across arrays
}
\description{
Plot trajectories of specific QC probes (e.g., biotin, cy3_hyb, housekeeping
gene probes, low stringency probes, etc.) across arrays
}
\usage{
plotQCCurves(
    esQC, 
    probes = c("biotin", "cy3_hyb", "housekeeping", 
      "low_stringency_hyb", "signal", "p95p05"), 
    labelVariable = "subjID",
    hybName = "Hybridization_Name",
    reporterGroupName = "Reporter_Group_Name",
    requireLog2 = TRUE, 
    projectName = "test", 
    plotOutPutFlag = FALSE, 
    cex = 1, 
    ylim = NULL, 
    xlab = "", 
    ylab = "intensity", 
    lwd = 3, 
    mar = c(10, 4, 4, 2) + 0.1,
    las = 2,
    cex.axis = 1,
    sortFlag = TRUE,
    varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat = c("\%m/\%d/\%Y", NA, NA),
    ...)
}
\arguments{
  \item{esQC}{
ExpressionSet object of QC probe profiles. \code{fData(esQC)} should
contains the variable \code{Reporter_Group_Name}. 
}
  \item{probes}{
A character vectors of QC probe names. By default, it includes
the following probe names \dQuote{biotin}, \dQuote{cy3_hyb}, 
  \dQuote{housekeeping}, \dQuote{low_stringency_hyb}, 
  \dQuote{signal}, \dQuote{p95p05}. For \dQuote{signal},
trajectories of 5th, 25th, 50th, 75th, and 95th percentiles of
the expression levels of all QC probes will be ploted. 
For \dQuote{p95p05}, the trajectory of the ratio of 95th percentile
to 5th percentile of the expression levels of all QC probes will be ploted.
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
\item{reporterGroupName}{
character string. indicating feature variable \code{Reporter_Group_Name} (QC probe's name).
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
}

  \item{\dots}{
Arguments to be passed to \code{\link{plot}}.
}
}
\value{
no return value.
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
    esQC.sim = genSimData.BayesNormal(nCpGs = 10, 
      nCases = 20, nControls = 20,
      mu.n = -2, mu.c = 2,
      d0 = 20, s02 = 0.64, s02.c = 1.5, testPara = "var",
      outlierFlag = FALSE, 
      eps = 1.0e-3, applier = lapply) 

    print(esQC.sim)

    fDat = fData(esQC.sim)
    esQC.sim$Hybridization_Name = sampleNames(esQC.sim)
    fDat$Reporter_Group_Name = c( rep("biotin", 5),
      rep("housekeeping", 5))
    fData(esQC.sim)=fDat

    # plot trajectories of the QC probes
    plotQCCurves(
      esQC = esQC.sim, 
      probes = c("biotin", "housekeeping"), 
      labelVariable = "subjID",
      hybName = "Hybridization_Name",
      reporterGroupName = "Reporter_Group_Name",
      requireLog2 = FALSE, 
      plotOutPutFlag = FALSE, 
      sortFlag = FALSE)

}

\keyword{ methods }
