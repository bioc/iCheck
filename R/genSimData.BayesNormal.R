# created on July 19, 2015
#   simulate data set with so that cases and controls
#     have different means and variances


# Phipson and Oshlack's (2014) simulation from null distribution
#
# for each CpG site, variance was first sampled from an
#   scaled inverse chi-squared distribution
#      sigma_i ~ scale-inv chi^2(d0, s0^2)
#  M value for each CpG was sampled from a normal distribution
#    with variance equal to the simulated variance
#

genSimData.BayesNormal=function(nCpGs, nCases, nControls,
  mu.n=-2, mu.c=2, d0=20, s02=0.64, s02.c=1.5, testPara="var", 
  outlierFlag = FALSE, eps=1.0e-3, applier=lapply) 
{
  tau2=s02
  v=d0
  # simulate variances for controls
  chi2=rchisq(n=nCpGs, df=v)
  sdVec = sqrt((v*tau2)/chi2)

  # simulate CpG from alternative hypothesis (un-equal variance)
  tau2.c=s02.c
  v=d0
  sdVec.c = sqrt((v*tau2.c)/chi2)

  m0=mu.n
  m1=mu.c

  ttLst=lapply(1:nCpGs, function(i) {
    x.c=rnorm(nCases, mean=mu.c, sd=sdVec.c[i])
    x.n=rnorm(nControls, mean=mu.n, sd=sdVec[i])

    v0=sdVec[i]
    v1=sdVec.c[i]
  
    if(testPara=="mean")
    {
       if(abs(m0-m1)<eps) 
       {
         memGenes=0
       } else {
         memGenes=1
       }
    } else if(testPara=="var") {
       if(abs(v0-v1)<eps) 
       {
         memGenes=0
       } else {
         memGenes=1
       }
    } else { # test for both mean difference and variance difference
       if(abs(m0-m1)<eps && abs(v0-v1)<eps) 
       {
         memGenes=0
       } else {
         memGenes=1
       }
    }
   
    return( c(memGenes, x.c, x.n) )
  })
  mat2=t(sapply(ttLst, function(x) {x}))
  memGenes=mat2[,1]
  mat=mat2[,-1]

  if(outlierFlag)
  {
    # add outlier
    datvec=c(mat)
    M.max=max(datvec, na.rm=TRUE)
  
    # randomly select one case to add outliers
    pos.outCase=sample(x=1:nCases, size=1, replace=FALSE)
    mat[, pos.outCase]=M.max
  }

  memSubj=c(rep(1, nCases), rep(0, nControls))

  nSubj=nCases+nControls
  subjID=paste("subj", 1:nSubj, sep="")
  probeID=paste("probe", 1:nCpGs, sep="")
  genes=paste("gene", 1:nCpGs, sep="")
  chr=rep(1, nCpGs)

  pDat=data.frame(
    arrayID=subjID,    
    memSubj=memSubj
  )
  rownames(pDat)=subjID

  fDat=data.frame(probe=probeID, gene=genes, chr=chr, memGenes=memGenes)
  rownames(fDat)=probeID

  rownames(mat)=probeID
  colnames(mat)=subjID

  # create ExpressionSet object
  aa<-new("AnnotatedDataFrame", data=pDat)
  #bb<-assayDataNew("lockedEnvironment", exprs=mat)
  #bb<-new("AssayData", exprs=mat)

  #es<-new("ExpressionSet",
  #    assayData=bb,
  #    phenoData = aa,
  #    annotation = "")

  es<-new("ExpressionSet",
      exprs=mat,
      phenoData = aa,
      annotation = "")

  Biobase::fData(es)=fDat

  invisible(es)

}

