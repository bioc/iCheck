#########
## create an ExpressionSet object
#########

genExprSet<-function(ex, pDat, fDat=NULL, 
  annotation="lumiHumanAll.db")
{
  # check the consistency of the column names of 'dat',
  # column names of 'fDat', and row names of 'pDat'
  cn.dat<-colnames(ex)
  rn.pdat<-rownames(pDat)

  aa<-match(cn.dat, rn.pdat)
  if(length(cn.dat)!=length(rn.pdat))
  {
    cat("Warning: No. of columns of ex=", length(cn.dat), "\n")
    cat("not equalt to that of pDat =", length(rn.pdat), "\n")
    diffxy<-setdiff(cn.dat, rn.pdat)
    if(length(diffxy))
    {
      cat("The sample in ex, but not in pDat are>>\n")
      print(diffxy)
      cat("\n")
    }
    diffyx<-setdiff(rn.pdat, cn.dat)
    if(length(diffyx))
    {
      cat("The sample in pDat, but not in ex are>>\n")
      print(diffyx)
      cat("\n")
    }
  }
  
  if(!any(is.na(aa)==TRUE))
  {
    pDat2<-pDat[aa,,drop=FALSE]
    identical(rownames(pDat2), colnames(ex))
    pDat3<-as(pDat2, "data.frame")
    aa<-new("AnnotatedDataFrame", data=pDat3)

    exprs<-as(ex, "matrix")
    es.raw<-new("ExpressionSet",
      exprs=exprs,
      phenoData = aa,
      annotation = annotation)
  } else {
    stop("Column names of ex != row names of pDat!\n")
  } 
  
  if(!is.null(fDat))
  { cn.fdat<-colnames(fDat)
    if(identical(sort(rownames(ex)), sort(rownames(fDat))))
    { 
      cc<-match(rownames(ex), rownames(fDat))
      Biobase::fData(es.raw)=fDat[cc,,drop=FALSE]
    } else {
      stop("Row names of ex != row names of dat.control!\n")
    }
  
  }
  
  invisible(es.raw)
}


