inputChecking<-function(pvalAdjMethod, es, probeID.var,
  gene.var, chr.var)
{
  pvalAdjMethod=match.arg(pvalAdjMethod, choices=c("holm", "hochberg", 
    "hommel", "bonferroni", "BH", "BY", "fdr", "none"))

  # extract feature data
  fDat<-fData(es)
  cn<-colnames(fDat)

  # get probe IDs, gene symbos and chromosome numbers
  pos1<-which(cn == probeID.var)
  if(length(pos1))
  { 
    probeIDs<-fDat[, pos1]
  } else {
    stop(paste("no variable", probeID.var, "in feature data", sep=" "))
  }

  pos1<-which(cn == gene.var)
  if(length(pos1))
  { 
    geneSymbols<-fDat[, pos1]
  } else {
    stop(paste("no variable", gene.var, "in feature data", sep=" "))
  }

  pos1<-which(cn == chr.var)
  if(length(pos1))
  { 
    chr<-fDat[, pos1]
  } else {
    stop(paste("no variable", chr.var, "in feature data", sep=" "))
  }


  # extract expression data
  dat<-exprs(es)

  # extract phenotype data
  pDat<-pData(es)

  res<-list(dat=dat, pDat=pDat, fDat=fDat,
    probeIDs = probeIDs, geneSymbols = geneSymbols, chr=chr)

}

getCovInterest<-function(formula, pos.var.interest, formulaType)
{
  if(pos.var.interest==0)
  {
    var.interest<-"(Intercept)"
  } else {
    fla<-as.character(formula)
    if(formulaType=="both")
    {
      # fla[1] = "~"
      # fla[2] = response variable
      # fla[3] = right-hand-side of "~" in formula
      fla3<-fla[3]
    } else {
      # fla[1] = "~"
      # fla[2] = right-hand-side of "~" in formula
      fla3<-fla[2]
    }

    tt<-strsplit(fla3, split="\\+")
    tt2<-as.character(sapply(tt, function(x) { limma::trimWhiteSpace(x)}))
    var.interest<-tt2[pos.var.interest]
  }
  return(var.interest)
}

# Check if the covariate of interest is a
# factor or interaction term with more than 2 levels.
# Return the columns in designMat containing all levels
#  of the covariate of interest.
#
checkCovariate<-function(
  formula, 
  cn.designMat, 
  pos.var.interest,
  formulaType="both")
{
  formulaType<-match.arg(formulaType, choices=c("both", "right"))

  var.interest<-getCovInterest(formula=formula, 
    pos.var.interest=pos.var.interest, formulaType=formulaType)

  if(pos.var.interest==0)
  { # the covariate of interest is the intercept
    # the intercept is the 1-st column of designMat
    colDesignMat<-1
    res<-list(var.interest=var.interest, colDesignMat=colDesignMat)
    return(res)
  } 

  cn<-cn.designMat

  # check if var.interest is an interaction term
  ttpos<-grep(pattern="\\*", x=var.interest, ignore.case =TRUE)
  ttpos2<-grep(pattern="\\:", x=var.interest, ignore.case =TRUE)
  if(length(ttpos) || length(ttpos2))
  { # var.interest is an interaction term
    if(length(ttpos))
    {
      tt3<-strsplit(var.interest, split="\\*")
    } else {
      tt3<-strsplit(var.interest, split="\\:")
    }
    # tt4 contains all covariates in the interaction term 
    tt4<-as.character(sapply(tt3, function(x) { limma::trimWhiteSpace(x)}))

    # check which columns of designMat contains all elements of tt4
    tt7<-unlist(sapply(cn, function(x) {
        tt5<-as.numeric(sapply(tt4, function(y) {
          tt6<-grep(pattern=y, x=x, fixed=TRUE)
          if(length(tt6))
          { 
            return(1) 
          } else {
            return(0) 
          }
        }))
        if(sum(tt5, na.rm=TRUE)==length(tt4))
        { # this column of designMat contains one level of var.interest  
          return(1)
        } else {
          # this column of designMat does not contain the leve of var.interest 
          return(0) 
        }
      }
    ))
    # return columns in designMat containing 
    # the all levels  of the covariate of the interest
    colDesignMat<-which(tt7==1)
    res<-list(var.interest=var.interest, colDesignMat=colDesignMat)
    return(res)
  } else {
    # not an interaction term
    # but might be a factor
    ttpos2<-grep(pattern=var.interest, x=cn, fixed=TRUE)
    ttlen<-length(ttpos2)

    if(ttlen==1)
    {
      # the covariate of interest is not a factor with
      # more than 2 levels
      # return columns in designMat containing 
      # the all levels  of the covariate of the interest
      colDesignMat<-ttpos2
      res<-list(var.interest=var.interest, colDesignMat=colDesignMat)
      return(res)
    } else {
      # each element of 'cn2' contains string var.interest 
      set<-1:length(cn)
      set2<-set[ttpos2]
      cn2<-cn[ttpos2]

      # check if cn2 contains interaction terms
      tt8<-as.numeric(unlist(sapply(cn2, function(x) {
          tt7<-grep(pattern="\\:", x=x, ignore.case=TRUE)
          if(length(tt7))
          {
            return(0) # not we want
          } else {
            return(1) # what we want
          }
        }
      )))
      # columns in designMat contains var.interest 
      ttpos3<-which(tt8 > 0)

      # return columns in designMat containing 
      # the all levels  of the covariate of the interest
      colDesignMat<-set2[ttpos3]
      res<-list(var.interest=var.interest, colDesignMat=colDesignMat)
      return(res)
    }
  }
}


getDesignMat.glm<-function(xi, pDat, formula)
{
  pDat.tmp<-pDat
  pDat.tmp$xi = xi
  designMat <- model.matrix(formula, pDat.tmp)
  invisible(designMat)
}

logit<-function(x)
{
  return(log(x)/(log(1-x)))
}


# put the variable interested in the first place of formula after ~
glmFitOneProbe<-function(
  xi, 
  pDat, 
  formula=xi~ratio, 
  logit=FALSE, 
  family=gaussian,
  verbose=FALSE)
{
  pDat.tmp<-pDat
  if(logit)
  {
    pDat.tmp$xi = log(xi/(1-xi))
  } else {
    pDat.tmp$xi = xi
  }

  # some times glm will exit abnormally 
  # even though we set na.action="na.exclude"
  # because some subjects containing missing values
  # Hence, we drop these subjects before calling glm
  pDat2<-pDat.tmp
  rownames(pDat2)<-1:nrow(pDat2)
  designMat <- model.matrix(formula, pDat2)
  nm<-colnames(designMat)
  rn<-as.numeric(rownames(designMat))
  # remove possible dropped subjects from data frame
  pDat.tmp<-pDat.tmp[rn,,drop=FALSE]

  ans<-try(summary(glm(formula=formula, family=family,
    data=pDat.tmp, na.action="na.exclude")))
  aaa<-attr(ans,which="class")
  if(length(aaa)>1)
  { aaa<-aaa[1] }
  if (aaa == "try-error")
  {
    nn<-length(nm)
    pvalVec <- rep(NA, nn)
    coefVec <- rep(NA, nn)
    statVec <- rep(NA, nn)
    cat("try-error; set pvalVec=NA and stat=NA\n")
  } else {

    if(verbose)
    { print(ans) }

    # some times, some covariates do not appear in
    # the ans$coefficients probably because
    # variances of these covariates are zero.
    # We need to pad NA for these covariates
    # to make sure the output for each gene probe 
    # has the same length
    nm2<-rownames(ans$coefficients)
    ttpos<-which(is.na(match(nm, nm2))==FALSE)
    statVec<-rep(NA, length(nm))
    statVec[ttpos]<-as.numeric(ans$coefficients[, 3])
    coefVec<-rep(NA, length(nm))
    coefVec[ttpos]<-as.numeric(ans$coefficients[, 1])
    pvalVec<-rep(NA, length(nm))
    pvalVec[ttpos]<-as.numeric(ans$coefficients[, 4])

  }

  res<-c(pvalVec, statVec, coefVec)
  names(res)<-c(paste("pval",nm, sep="."), paste("stat",nm, sep="."),
    paste("coef",nm, sep="."))
  return(res)
}



# pos.var.interest -- which covariate in the right-hand-side
# of ~ in the 'formula' is of the interest
#   pos.var.interest = 0 means intercept is of the interest
glmWrapper<-function(
  es,
  formula = FEV1~xi+age+gender, 
  pos.var.interest = 1,
  family=gaussian,
  logit=FALSE, 
  pvalAdjMethod="fdr", 
  alpha = 0.05, 
  probeID.var = "ProbeID", 
  gene.var = "Symbol", 
  chr.var = "Chromosome", 
  applier=lapply,
  verbose=TRUE 
  )
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

  t4<-applier(1:nrow(dat), function(i) {
      glmFitOneProbe(xi=dat[i,], pDat=pDat, formula=formula,
        logit=logit, family=family, verbose=FALSE)
    }
  )

  resMat<-t(sapply(t4, function(x) {x}))
  ttnc<-ncol(resMat)
  ttnc2<-ttnc/3

  # the first part is the p-value matrix
  pvalMat.unsorted<-resMat[,1:ttnc2, drop=FALSE]
  rownames(pvalMat.unsorted)<-featureNames(es)
  # the second part is the stat matrix
  statMat.unsorted<-resMat[,c((ttnc2+1):(2*ttnc2)), drop=FALSE]
  rownames(statMat.unsorted)<-featureNames(es)
  # the  third part is the coef matrix
  coefMat.unsorted<-resMat[,-c(1:(2*ttnc2)), drop=FALSE]
  rownames(coefMat.unsorted)<-featureNames(es)

  # check if the covariate of interest is a factor
  # or interaction term with more than 2 levels
  # first get example design matrix  
  designMat<-getDesignMat.glm(xi=dat[1,], pDat=pDat, 
    formula=formula)

  tmp<-checkCovariate(
    formula=formula, 
    cn.designMat=colnames(designMat), 
    pos.var.interest=pos.var.interest,
    formulaType="both")

  var.interest<-tmp$var.interest
  colDesignMat<-tmp$colDesignMat

  pval<-pvalMat.unsorted[, colDesignMat]
  stats<-statMat.unsorted[, colDesignMat]
  coef<-coefMat.unsorted[, colDesignMat]
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
    coef<-apply(coef, 1, function(x) { 
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

  res<-outputFunc2(
    pval=pval, 
    stats=stats, 
    coef=coef, 
    pvalAdjMethod=pvalAdjMethod,
    probeIDs=probeIDs, 
    geneSymbols=geneSymbols, 
    chr=chr, 
    pvalMat.unsorted=pvalMat.unsorted,
    statMat.unsorted=statMat.unsorted, 
    coefMat.unsorted=coefMat.unsorted, 
    dat=dat, 
    alpha=alpha,
    formula = formula,
    var.interest = var.interest,
    verbose=verbose)

  res$resMat=resMat

  invisible(res)
}

# pos.var.interest -- which covariate in the right-hand-side
# of ~ in the 'formula' is of the interest
#   pos.var.interest = 0 means intercept is of the interest
# likelihood ratio test
lkhrWrapper<-function(
  es,
  formulaReduced = FEV1~xi+gender, 
  formulaFull = FEV1~xi+age+gender, 
  family=gaussian,
  logit=FALSE, 
  pvalAdjMethod="fdr", 
  alpha = 0.05, 
  probeID.var = "ProbeID", 
  gene.var = "Symbol", 
  chr.var = "Chromosome", 
  applier=lapply,
  verbose=TRUE 
  )
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

  t4<-applier(1:nrow(dat), function(i) {
      lkhrFitOneProbe(xi=dat[i,], pDat=pDat, 
        formulaReduced=formulaReduced,
        formulaFull=formulaFull,
        logit=logit, family=family, verbose=FALSE)
    }
  )

  resMat<-t(sapply(t4, function(x) {x}))
  resFrame<-data.frame(
    probeIDs=probeIDs,
    geneSymbols=geneSymbols, 
    chr=chr, 
    Chisq=resMat[,1], 
    Df=resMat[,2], 
    pval=resMat[,3])
  resFrame$p.adj<-p.adjust(resFrame$pval, method="fdr")
  resFrame$pos<-1:nrow(resFrame)

  resFrame.s<-resFrame[order(resFrame$Chisq, decreasing = TRUE),]

  if(verbose)
  {
    cat("\n*********************************\n")
    nTop=min(c(20, nrow(resFrame.s)))
    cat("\nTop ", nTop, " tests>>>\n")
    print(resFrame.s[1:nTop, ])
    cat("\n")

    cat("\nformulaReduced>>\n")
    print(formulaReduced)
    cat("\nformulaFull>>\n")
    print(formulaFull)
    cat("\n")

    cat("\nNumber of tests>>>", nrow(resFrame.s), "\n")
    cat("\nNumber of arrays>>>", nrow(pDat), "\n")
    cat("\nNumber of tests with pvalue<0.05>>>", 
      sum(resFrame.s$pval<0.05, na.rm=TRUE), "\n")
    cat("\nNumber of tests with FDR adjusted pvalue<0.05>>>", 
      sum(resFrame.s$p.adj<0.05, na.rm=TRUE), "\n")

  }

  res<-list(frame=resFrame.s, frame.unsorted=resFrame)
  invisible(res)
}

# put the variable interested in the first place of formula after ~
lkhrFitOneProbe<-function(
  xi, 
  pDat, 
  formulaReduced=xi~ratio+gender, 
  formulaFull=xi~ratio+age+gender, 
  logit=FALSE, 
  family=gaussian,
  verbose=FALSE)
{
  pDat.tmp<-pDat
  if(logit)
  {
    pDat.tmp$xi = log(xi/(1-xi))
  } else {
    pDat.tmp$xi = xi
  }

  # some times glm will exit abnormally 
  # even though we set na.action="na.exclude"
  # because some subjects containing missing values
  # Hence, we drop these subjects before calling glm
  pDat2<-pDat.tmp
  rownames(pDat2)<-1:nrow(pDat2)
  designMatFull <- model.matrix(formulaFull, pDat2)
  nm<-colnames(designMatFull)
  rn<-as.numeric(rownames(designMatFull))
  # remove possible dropped subjects from data frame
  pDat.tmp<-pDat.tmp[rn,,drop=FALSE]

  ansReduced<-try((glm(formula=formulaReduced, family=family,
    data=pDat.tmp, na.action="na.exclude")))
  aaaReduced<-attr(ansReduced,which="class")

  if(length(aaaReduced)>1)
  { aaaReduced<-aaaReduced[1] }

  ansFull<-try((glm(formula=formulaFull, family=family,
    data=pDat.tmp, na.action="na.exclude")))
  aaaFull<-attr(ansFull,which="class")

  if(length(aaaFull)>1)
  { aaaFull<-aaaFull[1] }


  if (aaaReduced == "try-error" || aaaFull == "try-error")
  {
    cat("try-error; set pval=NA and stat=NA\n")
    res<-rep(NA, 3)
    names(res)<-c("Chisq", "Df", "pval")

  } else {

    if(verbose)
    { 
      cat("\nReduced model>>>\n")
      print(ansReduced) 

      cat("\nFull model>>>\n")
      print(ansFull) 
    }

    res.lr<-lrtest(ansReduced, ansFull)
    res<-c(res.lr$Chisq[2], res.lr$Df[2], res.lr$"Pr(>Chisq)"[2])
    names(res)<-c("Chisq", "Df", "pval")
    if(verbose)
    {
      cat("\nlikelihood ratio test>>>>\n")
      print(res.lr)
      cat("\n")
    }

  }

  #res<-c(pvalVec, statVec)
  #names(res)<-c(paste("pval",nm, sep="."), paste("stat",nm, sep="."))
  return(res)
}

outputFunc2<-function(pval, stats, coef, pvalAdjMethod,
  probeIDs, geneSymbols, chr, pvalMat.unsorted,
  statMat.unsorted, coefMat.unsorted, dat, alpha, formula, 
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
  mat<-cbind(pval, stats, coef, 1:nGenes)
  mat2<-mat[order(mat[,1]),]
  pos<-as.numeric(mat2[,4])
        
  # create output dataframe
  frame<-data.frame(probeIDs=probeIDs[pos], geneSymbols=geneSymbols[pos],
    chr=chr[pos],
    stats=mat2[,2], coef=mat2[,3], pval=mat2[,1], p.adj=p.adj[pos], pos=pos)
  rownames(frame)<-NULL

  # sorted p-value matrix
  pvalMat<-pvalMat.unsorted[pos, ,drop=FALSE]
  # sorted moderate t-test statistic matrix
  statMat<-statMat.unsorted[pos, ,drop=FALSE]
  # sorted coef matrix
  coefMat<-coefMat.unsorted[pos, ,drop=FALSE]

  frame.unsorted<-data.frame(probeIDs=probeIDs, geneSymbols=geneSymbols,
    chr=chr,
    stats=mat[,2], coef=mat[,3], pval=mat[,1], p.adj=p.adj, pos=1:nGenes)
  rownames(frame.unsorted)<-NULL
  rownames(statMat.unsorted)<-NULL
  rownames(pvalMat.unsorted)<-NULL
  rownames(coefMat.unsorted)<-NULL
  rownames(statMat)<-NULL
  rownames(coefMat)<-NULL
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
    coefMat.unsorted=coefMat.unsorted,
    pvalMat.unsorted=pvalMat.unsorted,
    memGenes=memGenes, memGenes2=memGenes2, mu1=mu1, mu2=mu2, mu3=mu3)

  invisible(res)

}

