
# make sure the first covariate in 'formula' is the variable of interest
# pos.var.interest indicates which covariate in right-hand-side
#   of "~" in 'formula' is the covariate of the interest
lmFitWrapper<-function(
  es, 
  formula=~as.factor(gender), 
  pos.var.interest = 1,
  pvalAdjMethod="fdr", 
  alpha=0.05, 
  probeID.var="ProbeID", 
  gene.var="Symbol", 
  chr.var="Chromosome", 
  verbose=TRUE)
{ 

  tt<-inputChecking(
    pvalAdjMethod=pvalAdjMethod, 
    es=es, 
    probeID.var=probeID.var,
    gene.var=gene.var, 
    chr.var=chr.var)

  dat<-tt$dat
  pDat<-tt$pDat
  fDat<-tt$fDat
  probeIDs=tt$probeIDs
  geneSymbols=tt$geneSymbols
  chr=tt$chr

  # some times designMat will loss some subjects
  # because these subjects containing missing values
  # Hence, we need to know which rows dropped
  pDat2<-pDat
  rownames(pDat2)<-1:nrow(pDat2)

  designMat <- model.matrix(formula, pDat2)
  rn<-as.numeric(rownames(designMat))
  # remove possible dropped subjects from expression data
  dat<-dat[,rn,drop=FALSE]
  rownames(designMat)<-colnames(dat)

  #####
  # gene differential analysis
  #####
  if(verbose)
  { 
    cat("dim(dat)>>\n"); print(dim(dat)); cat("\n");
    cat("Running lmFit...\n") 
  }
  fit = lmFit(dat, designMat)

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
  rownames(pvalMat.unsorted)<-featureNames(es)

  # unsorted moderate t-test statistic matrix
  statMat.unsorted <- ebFit$t
  if(is.null(dim(statMat.unsorted)))
  {
    statMat.unsorted<-matrix(statMat.unsorted, ncol=1)
    colnames(statMat.unsorted)<-"(Intercept)"
  }  

  rownames(statMat.unsorted)<-featureNames(es)
  

  # check if the covariate of interest is a factor
  # or interaction term with more than 2 levels
  tmp<-checkCovariate(
    formula=formula, 
    cn.designMat=colnames(designMat), 
    pos.var.interest=pos.var.interest,
    formulaType="right")

  var.interest<-tmp$var.interest
  colDesignMat<-tmp$colDesignMat

  pval<-pvalMat.unsorted[, colDesignMat]
  stats<-statMat.unsorted[, colDesignMat]
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
    dat=dat, 
    alpha=alpha,
    formula = formula,
    var.interest = var.interest,
    verbose=verbose)

  res$ebFit=ebFit 

  invisible(res)

}

outputFunc<-function(pval, stats, pvalAdjMethod,
  probeIDs, geneSymbols, chr, pvalMat.unsorted,
  statMat.unsorted, dat, alpha, formula, 
  var.interest,
  verbose)
{
  # quantiles of p-values for all covariates including intercept 
  pval.quantile<-apply(pvalMat.unsorted, 2, quantile, na.rm=TRUE)
  rownames(pval.quantile)<-c("min", "25%", "median", "75%", "max")
  colnames(pval.quantile)<-colnames(pvalMat.unsorted)

  # adjust p-value
  p.adj<-p.adjust(pval, method=pvalAdjMethod)

  # re-order genes by the ascending order of p-value
  # for the covariate of our interest
  nGenes<-length(pval)
  mat<-cbind(pval, stats, 1:nGenes)
  mat2<-mat[order(mat[,1]),]
  pos<-as.numeric(mat2[,3])
        
  # create output dataframe
  frame<-data.frame(probeIDs=probeIDs[pos], geneSymbols=geneSymbols[pos],
    chr=chr[pos],
    stats=mat2[,2], pval=mat2[,1], p.adj=p.adj[pos], pos=pos)
  rownames(frame)<-NULL

  # sorted p-value matrix
  pvalMat<-pvalMat.unsorted[pos, ,drop=FALSE]
  # sorted moderate t-test statistic matrix
  statMat<-statMat.unsorted[pos, ,drop=FALSE]

  frame.unsorted<-data.frame(probeIDs=probeIDs, geneSymbols=geneSymbols,
    chr=chr,
    stats=mat[,2], pval=mat[,1], p.adj=p.adj, pos=1:nGenes)
  rownames(frame.unsorted)<-NULL
  rownames(statMat.unsorted)<-NULL
  rownames(pvalMat.unsorted)<-NULL
  rownames(statMat)<-NULL
  rownames(pvalMat)<-NULL

  # number of significant gene probes after p-value adjustment
  n.sig<-sum(p.adj<alpha, na.rm=TRUE)
  if(verbose)
  { 
    nPrint=min(nrow(frame), 20)
    if(nPrint)
    { 
      print(frame[1:nPrint,])
    }

    cat("\npvalue quantiles for intercept and covariates>>\n")
    print(pval.quantile)
    cat("\n")

    cat("formula>>\n"); print(formula); cat("\n");
    cat("covariate of interest is ", var.interest, "\n")
    cat("Number of tests=", nrow(frame), "\n")
    cat("Number of arrays=", ncol(dat), "\n")
    cat("Number of significant tests (raw p-value < ", alpha, ")=", sum(frame$pval<0.05, na.rm=TRUE), "\n")
    cat("Number of significant tests after p-value adjustments=", n.sig, "\n")
    cat("\n\n**********************************************\n")
  }

  # the positive value of stats is over-expressed
  pos1<-which(p.adj<alpha & stats>0)
  len1<-length(pos1)

  # the negative value of stats is under-expressed
  pos3<-which(p.adj<alpha & stats<0)
  len3<-length(pos3)

  nGenes<-nrow(dat)
  memGenes<-rep(2, nGenes)
  memGenes2<-rep(0, nGenes)

  # mixing proportion
  if(len1)
  {
    pi.1<-len1/nGenes
    memGenes[pos1]<-1
  } else {
    pi.1<-0
  }
  if(len3)
  {
    pi.3<-length(pos3)/nGenes
    memGenes[pos3]<-3
  } else {
    pi.3<-0
  }
  pi.2<-1-pi.1-pi.3

  if(len1!=0 || len3!=0)
  {
    memGenes2[memGenes!=2]<-1
  } else {
    cat("No genes are differentially expressed!\n")
  }

  dat1<-dat[memGenes==1,,drop=FALSE]
  dat2<-dat[memGenes==2,,drop=FALSE]
  dat3<-dat[memGenes==3,,drop=FALSE]

  mu1<-apply(dat1, 2, mean, na.rm=TRUE)
  mu2<-apply(dat2, 2, mean, na.rm=TRUE)
  mu3<-apply(dat3, 2, mean, na.rm=TRUE)

  res<-list(n.sig=n.sig, frame=frame, 
    pvalMat=pvalMat,
    pval.quantile=pval.quantile,
    frame.unsorted=frame.unsorted,
    statMat.unsorted=statMat.unsorted,
    pvalMat.unsorted=pvalMat.unsorted,
    memGenes=memGenes, memGenes2=memGenes2, mu1=mu1, mu2=mu2, mu3=mu3)

  invisible(res)

}

