# intersection strategy
myCombine<-function(x, y, combineRule="intersection", featureID.es="Probe_Sequence", suffix="_2")
{
  combineRule = match.arg(combineRule, choices=c("union", "intersection"))
  if(combineRule=="intersection")
  {
    res<-myCombineIntersect(x, y, featureID.es=featureID.es, suffix=suffix)
  } else {
    res<-myCombineUnion(x, y, featureID.es=featureID.es, suffix=suffix)
  }
  return(res)
}

myCombineIntersect<-function(x, y, featureID.es="Probe_Sequence", suffix="_2")
{
        if(!identical(featureNames(x), as.character(fData(x)[, featureID.es])))
        {
          cat("Error: featureNames(x) not equal to fData(x)[, featureID.es]!\n")
          return(NULL)
        } 
        if (missing(x)) {
          cat("Error: x is missing!\n")
          return(NULL)
        }
        if (missing(y)) {
          cat("Error: y is missing!\n")
          return(NULL)
        }
        if (class(x) != class(y))
        {
                cat(paste("objects must be the same class, but are ",
                 class(x), ", ", class(y), sep=""))
                 return(NULL)
        }

        if (!identical(sort(featureNames(x)), sort(featureNames(y)))) 
        {
          cat('Two objects have different feature names.\n')
        }
        ## determine whether there are duplicated sample names
        sampleName.x <- sampleNames(x)
        sampleName.y <- sampleNames(y)
        int<-intersect(sampleName.x, sampleName.y)
        if (length(int)) {
                cat('Error: Two objects have some common sample names!\n')
                return(NULL)

        }

        # merge duplicated sequences
        x<-mergeDuplicates(x, featureID.es)
        y<-mergeDuplicates(y, featureID.es)

        featureName.x <- featureNames(x)
        featureName.y <- featureNames(y)
        featureName.com <- intersect(featureName.x, featureName.y)
        len.com<-length(featureName.com)
        if (len.com) {
                if(len.com < length(featureName.x) || len.com < length(featureName.y))
                {
                  warning('Two objects have different featureNames, only the common ones were used.')
                }
        } else {
                cat('Error: Two objects have totally different featureNames!\n')
                return(NULL)
        }
        #x <- x[which(featureName.x %in% featureName.com),]
        x<-x[match(featureName.com, featureName.x),]
        #y <- y[which(featureName.y %in% featureName.com),]
        y<-y[match(featureName.com, featureName.y),]
        
        ## make sure two objects have the sample order of features
        if(!identical(sort(featureNames(x)), sort(featureNames(y))))
        {
          cat('Error: sort(featureNames(x)) not equal to sort(featureNames(y))!\n')
          return(NULL)
        }
        pos<-match(featureNames(x), featureNames(y))
        y<-y[pos,]

        history.submitted <- as.character(Sys.time())
        dimm.x <- dim(x) 
        dimm.y <- dim(y)

        x.bak<-x
        # combine x and y and update x
        Biobase::phenoData(x) <- Biobase::combine(phenoData(x.bak),phenoData(y))
        Biobase::assayData(x) <- Biobase::combine(Biobase::assayData(x.bak), Biobase::assayData(y))

        # make sure the sampleNames(x) equal to rownames of phenoData(x)
        sn<-sampleNames(x)
        rn<-rownames(pData(x))
        if(!identical(sort(sn), sort(rn)))
        {
          cat("Error: sort(sampleNames(x)) not equal to sort(rownames(pData(x)))!\n")
          return(NULL)
        }
        if(!identical(sn, rn))
        {
          pos<-match(sn, rn)
          pDat<-pData(x)
          pDat<-pDat[pos,]
          rn<-rn[pos]
          rownames(pDat)<-rn
          Biobase::pData(x)<-pDat
        }

        chipver.x<-unique(as.character(x$Chip_Manifest))
        chipver.y<-unique(as.character(y$Chip_Manifest))
        if(!is.null(chipver.x) && !is.null(chipver.y) && length(chipver.x) > 0 && length(chipver.y) > 0 && !identical(chipver.x, chipver.y))
        {
          # combine feature data
          cat("Warning: chipver.x not equal to chipver.y!\n")
          fDat.x<-fData(x)
          fDat.y<-fData(y)
          cn.x<-colnames(fDat.x)
          cn.y<-colnames(fDat.y)
          pos.y<-which(colnames(fData(y))==featureID.es)
          cn.y<-paste(cn.y, suffix, sep="")
          colnames(fDat.y)<-cn.y
          colnames(fDat.y)[pos.y]<-featureID.es
          Biobase::fData(y)<-fDat.y
          #fData(x)$mysid<-as.character(fData(x)[, featureID.es])
          #fData(y)$mysid<-as.character(fData(y)[, featureID.es])
          Biobase::featureData(x) <- Biobase::combine(Biobase::featureData(x.bak), Biobase::featureData(y))

          # make sure the featureNames(x) equal to rownames of featureData(x)
          sn<-featureNames(x)
          rn<-rownames(fData(x))
          if(!identical(sort(sn), sort(rn)))
          {
            cat("Error: sort(featureNames(x)) not equal to sort(rownames(fData(x)))!\n")
            return(NULL)
          }
          if(!identical(sn, rn))
          {
            pos<-match(sn, rn)
            fDat<-fData(x)
            fDat<-fDat[pos,]
            rn<-rn[pos]
            rownames(fDat)<-rn
            Biobase::fData(x)<-fDat
          }
        } 
    
        Biobase::experimentData(x) <- Biobase::combine(Biobase::experimentData(x.bak), Biobase::experimentData(y))
        Biobase::protocolData(x) <- Biobase::combine(Biobase::protocolData(x.bak),Biobase::protocolData(y))

        ## combining the QC information
        if (length(x@QC) > 0 && length(y@QC) > 0) {
                if (!is.null(x@QC$BeadStudioSummary) && !is.null(y@QC$BeadStudioSummary)) {
                        if (ncol(x@QC$BeadStudioSummary) == ncol(y@QC$BeadStudioSummary) && ncol(x@QC$BeadStudioSummary) > 0) {
                                BeadStudioSummary <- rbind(x@QC$BeadStudioSummary, y@QC$BeadStudioSummary)
                                x@QC$BeadStudioSummary <- BeadStudioSummary
                        } else {
                                x@QC <- list()
                        }
                } else {
                        x@QC <- list()
                }
                if (!is.null(x@QC$sampleSummary) && !is.null(y@QC$sampleSummary)) {
                        if (nrow(x@QC$sampleSummary) == nrow(y@QC$sampleSummary) && nrow(x@QC$sampleSummary)>0) {
                                sampleSummary <- cbind(x@QC$sampleSummary, y@QC$sampleSummary)
                                x@QC$sampleSummary <- sampleSummary
                        } else {
                                x@QC <- list()
                        }
                } else {
                        x@QC <- list()
                }
                if (!is.null(x@QC)) {
                        history.x <- x@QC$history
                        if (is.null(history.x)) history.x <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
                        if (is.null(history.x$lumiVersion)) history.x$lumiVersion <- rep(NA, nrow(history.x))
                        history.y <- y@QC$history
                        if (is.null(history.y)) history.y <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
                        if (is.null(history.y$lumiVersion)) history.y$lumiVersion <- rep(NA, nrow(history.y))
                        x@QC$history <- rbind(history.x, history.y)
                }
        } else {
                x@QC <- list()
        }

        ## VST transformation parameters
        if (!is.null(attr(x, 'vstParameter')) && !is.null(attr(y, 'vstParameter'))) {
                vstParameter.x <- attr(x, 'vstParameter')
                vstParameter.y <- attr(y, 'vstParameter')
                if (is.null(nrow(vstParameter.x))) {
                        vstParameter.x <- matrix(vstParameter.x, nrow=1)
                }
                if (is.null(nrow(vstParameter.y))) {
                        vstParameter.y <- matrix(vstParameter.y, nrow=1)
                }
                if (nrow(vstParameter.x) != dimm.x[2] || nrow(vstParameter.y) != dimm.y[2]) {
                        attr(x, 'vstParameter') <- attr(x, 'transformFun') <- NULL
                } else {
                        attr(x, 'vstParameter') <- rbind(attr(x, 'vstParameter'), attr(y, 'vstParameter'))
                        attr(x, 'transformFun') <- c(attr(x, 'transformFun'), attr(y, 'transformFun'))
                }
        }

        ## controlData information
        if (nrow(x@controlData) > 0 && nrow(y@controlData) > 0) {
          control.x<-x@controlData
          control.y<-y@controlData

          ## determine whether there are duplicated sample names
          sampleName.x <- colnames(control.x)
          sampleName.y <- colnames(control.y)
          int<-intersect(sampleName.x, sampleName.y)
          if (length(int)) {
                  cat("Error: Two controlData objects have some duplicated sample names!\n") 
                  return(NULL)
          }
          featureName.x <- rownames(control.x)
          featureName.y <- rownames(control.y)
          featureName.com <- intersect(featureName.x, featureName.y)
          len.com<-length(featureName.com)
          if (len.com > 0) {
             warning('Two controlData objects have different rownames, only the common ones were used.')
             #control.x <- control.x[which(rownames(control.x) %in% featureName.com),] 
             control.x <- control.x[match(featureName.com, rownames(control.x)),] 
             #control.y <- control.y[which(rownames(control.y) %in% featureName.com),]
             control.y <- control.y[match(featureName.com, rownames(control.y)),]
             ## make sure two objects have the sample order of features
             if(identical(sort(rownames(control.x)), sort(rownames(control.y))))
             {
                pos<-match(rownames(control.x), rownames(control.y))
                control.y<-control.y[pos,]
     
                control.x <- Biobase::combine(control.x, control.y)
                if(!identical(colnames(control.x), sampleNames(x)))
                {
                  cat("Error: colnames(control.x) not equal to sampleNames(x)!\n")
                  return(NULL)
                }
                x@controlData<-control.x

             } else {
               cat('Warning: sort(rownames(control.x)) not equal to sort(rownames(control.y))!\n')
               #x@controlData<-NULL
               x@controlData<-data.frame()
               return(NULL)
             }

          } else {
             cat('Warning: Two controlData objects have totally different rownames!\n')
             #x@controlData<-NULL
             x@controlData<-data.frame()
          }
        } else {
          #x@controlData <- NULL
          x@controlData<-data.frame()
        }

        # history tracking
        history.finished <- as.character(Sys.time())
        history.command <- capture.output(print(match.call(combine)))
        x@history<- rbind(x.bak@history, y@history)
        if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
        lumiVersion <- packageDescription('lumi')$Version
        x@history<- rbind(x@history,
               data.frame(submitted=history.submitted, finished=history.finished, command=history.command, lumiVersion=lumiVersion))
        #featureNames(x)<-as.character(fData(x)$Probe_Sequence)
        return(x)
}

# union strategy
myCombineUnion<-function(x, y, featureID.es="Probe_Sequence", suffix="_2")
{
  if(!identical(featureNames(x), as.character(fData(x)[, featureID.es])))
  {
    cat("Error: featureNames(x) not equal to fData(x)[, featureID.es]!\n")
    return(NULL)
  } 
  if (missing(x)) {
    cat("Error: x is missing!\n")
    if(!missing(y))
    {
      return(y) 
    } 
    return(NULL)
  }
  if (missing(y)) {
    if(!missing(x))
    {
      return(x)
    } else {
      return(NULL)
    }
  }

  if (class(x) != class(y))
  {
          cat(paste("objects must be the same class, but are ",
           class(x), ", ", class(y), sep=""))
          return(NULL)
  }

  if (!identical(sort(featureNames(x)), sort(featureNames(y)))) 
  {
    cat('Warning: Two data sets have different row names.\n')
  }
  ## determine whether there are duplicated sample names
  sampleName.x <- sampleNames(x)
  sampleName.y <- sampleNames(y)

  int<-intersect(sampleName.x, sampleName.y)
  if (length(int)) {
          cat('Error: Two objects have some duplicated sample names!\n')
          return(NULL)
  }

  # merge duplicated sequences
  x<-mergeDuplicates(x, featureID.es)
  seq.x <- as.character(fData(x)[, featureID.es])

  y<-mergeDuplicates(y, featureID.es)
  seq.y <- as.character(fData(y)[, featureID.es])

  seq.com<-intersect(seq.x, seq.y)
  if(length(seq.com))
  {
    if(!identical(sort(seq.x), sort(seq.y)))
    {
      warning('Two data sets have different probe sequence. The union were used.')
    }
  } else {
    cat('Two data sets have totally different featureNames. The union were used\n')
  }

  history.submitted <- as.character(Sys.time())
  dimm.x <- dim(x) 
  dimm.y <- dim(y)


  ## combine pheno data
  phenoData.xy <- Biobase::combine(phenoData(x),phenoData(y))

  # merge feature data
  chipver.x<-unique(as.character(x$Chip_Manifest))
  chipver.y<-unique(as.character(y$Chip_Manifest))
  if (is.null(chipver.y) || is.null(chipver.x) || length(chipver.x) == 0 || length(chipver.y) == 0 || chipver.y %in% chipver.x)  
  {
    fDat.xy<-fData(x)
  } else { 
    # combine feature data
    cat("Warning: chipver.x not equal to chipver.y!\n")
    fDat.x<-fData(x)
    fDat.y<-fData(y)
    cn.x<-colnames(fDat.x)
    cn.y<-colnames(fDat.y)
    pos.y<-which(colnames(fData(y))==featureID.es)
    cn.y<-paste(cn.y, suffix, sep="")
    colnames(fDat.y)<-cn.y
    colnames(fDat.y)[pos.y]<-featureID.es
    Biobase::fData(y)<-fDat.y
    fDat.xy<-merge(x=fData(x), y=fData(y), by=featureID.es,
      all=TRUE, sort=FALSE, suffixes = c("",suffix))
  }

  eltVec.x<-Biobase::assayDataElementNames(x)
  eltVec.y<-Biobase::assayDataElementNames(y)

  elt.xy<-intersect(eltVec.x, eltVec.y)
  len.elt.xy<-length(elt.xy)
  if(len.elt.xy==0)
  {
    cat("No overlapping arrayDataElementNames between 'x' and 'y'!\n")
    return(NULL)
  }

  elt.xnoty<-setdiff(eltVec.x, eltVec.y)
  if(length(elt.xnoty))
  {
    cat("\nThe following assayDataElementNames are in 'x', but not in 'y'!\n")
    print(elt.xnoty)
    cat("\nThese elements will not be in the combined object!\n") 
  }

  elt.ynotx<-setdiff(eltVec.y, eltVec.x)
  if(length(elt.ynotx))
  {
    cat("\nThe following assayDataElementNames are in 'y', but not in 'x'!\n")
    print(elt.ynotx)
    cat("\nThese elements will not be in the combined object!\n") 
  }
  storage.mode.x <- Biobase::storageMode(x)
  if(storage.mode.x == "lockedEnvironment")
  {
    aData.x<-Biobase::copyEnv(x@assayData)
  } else {
    aData.x <- x@assayData
  }

  storage.mode.y <- Biobase::storageMode(y)
  if(storage.mode.y == "lockedEnvironment")
  {
    aData.y<-Biobase::copyEnv(y@assayData)
  } else {
    aData.y <- y@assayData
  }

  myassayData<-list()
  if(length(which(elt.xy %in% c("exprs", "se.exprs")))==2)
  {
    cmd <- 'xy <- new("LumiBatch"'
  } else if (length(which(elt.xy %in% c("exprs")))==1) {
    cmd <- 'xy <- new("ExpressionSet"'
  } else {
    cat("elt.xy does not contain 'exprs'!\n")
    return(NULL)
  }
 
  # merge each common element
  for(elt in elt.xy)
  {
    aD.x.elt<-as.data.frame(aData.x[[elt]])
    aD.y.elt<-as.data.frame(aData.y[[elt]])
    
    int.xy<-intersect(colnames(aD.x.elt), colnames(aD.y.elt))
    len.int.xy<-length(int.xy)
    if(len.int.xy)
    {
      cat("Error: the following samples are in both x and y>>\n");
      print(int.xy);
      cat("\n");
      return(NULL)
    }

    aD.x.elt$myid<-as.character(fData(x)[, featureID.es])
    aD.y.elt$myid<-as.character(fData(y)[, featureID.es])
  
    aD.xy.elt<-merge(x=aD.x.elt, y=aD.y.elt, by="myid", all=TRUE, sort=FALSE)

    # make sure the order of rows is the same as that of fDat.xy
    aa<-match(as.character(fDat.xy[, featureID.es]), as.character(aD.xy.elt$myid))
    if(any(is.na(aa)==TRUE))
    {
      cat("not all", featureID.es, " in fDat.xy are also in", elt, "\n")
      return(NULL)
    }
    aD.xy.elt<-aD.xy.elt[aa,]

    rn<-as.character(aD.xy.elt$myid)
    rownames(aD.xy.elt)<-rn

    if(!identical(as.character(fDat.xy[, featureID.es]), rn))
    {
      cat("After combining, aD.xy.elt$myid != fData(xy)[, ", featureID.es, "]!\n")
      return(NULL)
    }

    ttpos<-which(colnames(aD.xy.elt) == "myid")
    aD.xy.elt<-aD.xy.elt[,-ttpos, drop=FALSE]

    myassayData[[elt]]<-as.matrix(aD.xy.elt)
    cmd <- paste(cmd, ',', elt, '=myassayData[["', elt, '"]]', sep="")
  }

  ### combine pheno data
  #phenoData.xy <- Biobase::combine(phenoData(x),phenoData(y))

  pos<-match(colnames(myassayData[["exprs"]]), sampleNames(phenoData.xy))
  if(any(is.na(pos)==TRUE))
  {
    cat("Error: sampleNames(xy) != sampleNames(phenoData.xy)!\n")
    return(NULL)
  }
  if(dim(phenoData.xy)[2]==1)
  {
    phenoData.xy2<-phenoData.xy[pos]
  } else {
    phenoData.xy2<-phenoData.xy[pos,]
  }

  cmd <- paste(cmd, ', phenoData=phenoData.xy2)', sep="")

  #cmd <- paste(cmd, ')')
  eval(parse(text=cmd))

  ## add feature data
  rownames(fDat.xy)<-as.character(fDat.xy[, featureID.es])
  bb<-new("AnnotatedDataFrame", data=fDat.xy)
  Biobase::featureData(xy)<-bb
  #fData(xy)<-fDat.xy

  ## combine experimentData
  Biobase::experimentData(xy) <- Biobase::combine(Biobase::experimentData(x),Biobase::experimentData(y))

  ## combine protocolData
  Biobase::protocolData(xy) <- Biobase::combine(Biobase::protocolData(x),Biobase::protocolData(y))

  ## combining the QC information
  if (length(x@QC) > 0 && length(y@QC) > 0) {
          if (!is.null(x@QC$BeadStudioSummary) && !is.null(y@QC$BeadStudioSummary)) {
                  if (ncol(x@QC$BeadStudioSummary) == ncol(y@QC$BeadStudioSummary) && ncol(x@QC$BeadStudioSummary) > 0) {
                          BeadStudioSummary <- rbind(x@QC$BeadStudioSummary, y@QC$BeadStudioSummary)
                          xy@QC$BeadStudioSummary <- BeadStudioSummary
                  } else {
                          xy@QC <- list()
                  }
          } else {
                  xy@QC <- list()
          }
          if (!is.null(x@QC$sampleSummary) && !is.null(y@QC$sampleSummary)) {
                  if (nrow(x@QC$sampleSummary) == nrow(y@QC$sampleSummary)) {
                          sampleSummary <- cbind(x@QC$sampleSummary, y@QC$sampleSummary)
                          xy@QC$sampleSummary <- sampleSummary
                  } else {
                          xy@QC <- list()
                  }
          } else {
                  xy@QC <- list()
          }
          if (!is.null(x@QC)) {
                  history.x <- x@QC$history
                  if (is.null(history.x)) history.x <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
                  if (is.null(history.x$lumiVersion)) history.x$lumiVersion <- rep(NA, nrow(history.x))
                  history.y <- y@QC$history
                  if (is.null(history.y)) history.y <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
                  if (is.null(history.y$lumiVersion)) history.y$lumiVersion <- rep(NA, nrow(history.y))
                  xy@QC$history <- rbind(history.x, history.y)
          }
  } else {
          xy@QC <- list()
  }

  ## VST transformation parameters
  if (!is.null(attr(x, 'vstParameter')) && !is.null(attr(y, 'vstParameter'))) {
          vstParameter.x <- attr(x, 'vstParameter')
          vstParameter.y <- attr(y, 'vstParameter')
          if (is.null(nrow(vstParameter.x))) {
                  vstParameter.x <- matrix(vstParameter.x, nrow=1)
          }
          if (is.null(nrow(vstParameter.y))) {
                  vstParameter.y <- matrix(vstParameter.y, nrow=1)
          }
          if (nrow(vstParameter.x) != dimm.x[2] || nrow(vstParameter.y) != dimm.y[2]) {
                  attr(xy, 'vstParameter') <- attr(x, 'transformFun') <- NULL
          } else {
                  attr(xy, 'vstParameter') <- rbind(attr(x, 'vstParameter'), attr(y, 'vstParameter'))
                  attr(xy, 'transformFun') <- c(attr(x, 'transformFun'), attr(y, 'transformFun'))
          }
  }

  ## controlData information
  if (nrow(x@controlData) > 0 && nrow(y@controlData) > 0) {
    control.x<-x@controlData
    control.y<-y@controlData

    ## determine whether there are duplicated sample names
    sampleName.x <- colnames(control.x)
    sampleName.y <- colnames(control.y)
    int<-intersect(sampleName.x, sampleName.y)
    if (length(int)) {
            cat("Error: Two controlData objects have some duplicated sample names!\n") 
            return(NULL)
    }

    if(identical(rownames(control.x), rownames(control.y)))
    { 
      rn<-rownames(control.x)
      tt<-cbind(as.matrix(control.x), as.matrix(control.y)) 
      rownames(tt)<-rn
      xy@controlData<-as.data.frame(tt) 

      pos<-match(sampleNames(xy), colnames(xy@controlData))
      if(any(is.na(pos)==TRUE))
      {
        cat("Error: sampleNames(xy) != colnames(xy@controlData)!\n")
        return(NULL)
      }
      xy@controlData<-xy@controlData[,pos]

    } else if (length(intersect(colnames(control.x), colnames(control.y)))) {
      cat("Error: control.x and control.y have common samples!\n")
      return(NULL)

    } else {
      rn.x<-row.names(control.x)
      rn.y<-row.names(control.y)
      common<-intersect(rn.x, rn.y)
      if(length(common))
      {
        #control.x2<-control.x[which(rn.x %in% common),]
        control.x2<-control.x[match(common, rn.x),]
        #control.y2<-control.y[which(rn.y %in% common),]
        control.y2<-control.y[match(common, rn.y),]
        xy@controlData<-as.data.frame(cbind(control.x2, control.y2)) 

        pos<-match(sampleNames(xy), colnames(xy@controlData))
        if(any(is.na(pos)==TRUE))
        {
          cat("Error: sampleNames(xy) != colnames(xy@controlData)!\n")
          return(NULL)
        }
        xy@controlData<-xy@controlData[,pos]

      } else { 
        #xy@controlData <- NULL
        xy@controlData <- data.frame()
      }
    }
  } else {
    #xy@controlData <- NULL
    xy@controlData <- data.frame()
  }

  # history tracking
  history.finished <- as.character(Sys.time())
  history.command <- capture.output(print(match.call(combine)))
  xy@history<- rbind(x@history, y@history)
  if (is.null(xy@history$lumiVersion)) xy@history$lumiVersion <- rep(NA, nrow(x@history))
  lumiVersion <- packageDescription('lumi')$Version
  history<- rbind(xy@history,
         data.frame(submitted=history.submitted, finished=history.finished, command=history.command, lumiVersion=lumiVersion))
  xy@history<-history
  
  Biobase::featureNames(xy)<-as.character(fData(xy)$Probe_Sequence)
  return(xy)
}


