# gene-wise linear model for log2 differences
# By default, formula = ~1
# The intercept is what we are interested in.
# The non-nominal covariates should be mean-centered
# to get correct effect size for the intercept.

lmFitPaired<-function(
  esDiff,
  formula=~1,
  pos.var.interest = 0,
  pvalAdjMethod="fdr", 
  alpha=0.05, 
  probeID.var="ProbeID", 
  gene.var="Symbol", 
  chr.var="Chromosome",
  verbose = TRUE)
{
  tt<-inputChecking(
    pvalAdjMethod=pvalAdjMethod, 
    es=esDiff, 
    probeID.var=probeID.var,
    gene.var=gene.var, 
    chr.var=chr.var)

  diffDat<-tt$dat
  pDat<-tt$pDat
  fDat<-tt$fDat
  probeIDs=tt$probeIDs
  geneSymbols=tt$geneSymbols
  chr=tt$chr

  # some times designMat will lose some subjects
  # because these subjects containing missing values
  # Hence, we need to know which rows dropped
  pDat2<-pDat
  rownames(pDat2)<-1:nrow(pDat2)

  if(is.null(formula))
  {
    designMat<-matrix(rep(1, ncol(diffDat)), ncol=1)
    rownames(designMat)<-1:nrow(designMat)
    colnames(designMat)<-"(Intercept)"
  } else {
    flist = as.list(formula)
    if (length(flist) > 1)
    {
      if(flist[[2]]!=1)
      { # get design matrix
        designMat<-model.matrix(formula, data=pDat2)
      } else {
        designMat<-matrix(rep(1, ncol(diffDat)), ncol=1)
        rownames(designMat)<-1:nrow(designMat)
        colnames(designMat)<-"(Intercept)"
      }
    } else {
      cat("Examples of formula are ~1 or ~cov1+cov2\n")
      stop("formula not correct.\n")
    }
  }

  rn<-as.numeric(rownames(designMat))
  # remove possible dropped subjects from expression data
  diffDat<-diffDat[,rn,drop=FALSE]
  rownames(designMat)<-colnames(diffDat)

  if(verbose)
  {
    cat("Running lmFit for paired data ...\n")
  }
  fit = lmFit(diffDat, designMat)

  if(verbose)
  { cat("Running eBayes...\n") }
  ebFit = eBayes(fit)

  ####
  # prepare output
  ####
  if(verbose)
  { cat("Preparing output...\n") }

  # unsorted pvalue matrix
  pvalMat.unsorted<-ebFit$p.value
  if(is.null(dim(pvalMat.unsorted)))
  {
    pvalMat.unsorted<-matrix(pvalMat.unsorted, ncol=1)
    colnames(pvalMat.unsorted)<-"(Intercept)"
  }  
  rownames(pvalMat.unsorted)<-featureNames(esDiff)

  # unsorted moderate t-test statistic matrix
  statMat.unsorted <- ebFit$t
  if(is.null(dim(statMat.unsorted)))
  {
    statMat.unsorted<-matrix(statMat.unsorted, ncol=1)
    colnames(statMat.unsorted)<-"(Intercept)"
  }  

  rownames(statMat.unsorted)<-featureNames(esDiff)

  # check if the covariate of interest is a factor
  # or interaction term with more than 2 levels
  tmp<-checkCovariate(
    formula=formula, 
    cn.designMat=colnames(designMat), 
    pos.var.interest=pos.var.interest,
    formulaType="right")

  var.interest<-tmp$var.interest
  colDesignMat<-tmp$colDesignMat

  if(is.null(dim(pvalMat.unsorted)))
  {
    pvalMat.unsorted<-matrix(pvalMat.unsorted, ncol=1)
    statMat.unsorted<-matrix(statMat.unsorted, ncol=1)
  }  
  pval<-pvalMat.unsorted[, colDesignMat, drop=FALSE]
  stats<-statMat.unsorted[, colDesignMat, drop=FALSE]

  len<-length(colDesignMat)
  if(len>1)
  {
    # take the minimum p-value and its test statistic 
    pval<-apply(pval, 1, min, na.rm=TRUE)
    stats<-apply(stats, 1, function(x) { 
      absx<-abs(x)
      pos<-which(absx == max(absx, na.rm=TRUE))
      return(x[pos[1]])
    })
    if(verbose)
    {
      cat("\n************************\n")
      cat("Warning: the covariate of interest ", var.interest, " is a factor\n")
      cat("or interaction term with more than 2 levels!\n") 
      cat("Only the smallest p-value will be output!\n")
      cat("Better to run the function lkhWrapper!\n")
      cat("************************\n")
    }
  }

  res<-outputFunc(
    pval=pval, 
    stats=stats, 
    pvalAdjMethod=pvalAdjMethod,
    probeIDs=probeIDs, 
    geneSymbols=geneSymbols, 
    chr=chr, 
    pvalMat.unsorted=pvalMat.unsorted,
    statMat.unsorted=statMat.unsorted, 
    dat=diffDat, 
    alpha=alpha,
    formula = formula,
    var.interest = var.interest,
    verbose=verbose)

  res$ebFit=ebFit 

  invisible(res)
}

