###########
# modified on Nov. 17, 2016
#  (1) revised 'scatterPlots', 'boxPlots' so that
#      the user can choose not output 'p.adj' in the subtitle
#  (2) add input options 'outcomeFlag' and 'fitLineFlag' to 'scatterPlots'
#
# modified on March 31, 2016
#  (1) added functions 'scatterPlots', 'boxPlots', 'densityPlots'
#
############
# modified on Aug 20, 2015
#  (1) deleted un-used functions: pcaPlotFunc
#
# modified on Feb 14, 2014
#  (1) fixed a bug in  'getPCAFunc': 
#       forgot updating 'es' after deleting probes with missing values
#  
# modified on Jan 11, 2014
#  (1) revised 'getPCAFunc': set 'corFlag = FALSE'
#
# modified on Sept 29, 2013
#  (1) revised 'getPCAFunc': added line
#      dat=na.omit(dat)
#    right before the calling of the function prcompt

# modified on July 21, 2010
#  (1) added argument 'labelVariable' to the functions 'plotQCCurves'
#       and 'plotCurves'
# created on July 3, 2010
#  (1) read in data
#
###########

# plot trajectories of QC probe expression levels across arrays
# for specific QC probes: biotin, cy3_hyb, housekeeping, 
# low_stringency_hyb, etc. 
plotQCCurves<-function(esQC, 
    probes=c("biotin", "cy3_hyb", "housekeeping",
    "low_stringency_hyb", "signal", "p95p05"), 
    labelVariable="subjID",
    hybName = "Hybridization_Name",
    reporterGroupName = "Reporter_Group_Name",
    requireLog2=TRUE, 
    projectName="test", 
    plotOutPutFlag=FALSE, 
    cex=1, 
    ylim=NULL, 
    xlab="", 
    ylab="intensity", 
    lwd=3, 
    mar=c(10, 4, 4, 2) + 0.1,
    las=2,
    cex.axis=1,
    sortFlag = TRUE,
    varSort=c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat=c("%m/%d/%Y", NA, NA),
    ...)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }
  cat("probes>>\n"); print(probes); cat("\n");

  if(sortFlag)
  { esQC<-sortExpressionSet(esQC, varSort, timeFormat) }

  datQC<-getDatColnames(esQC, hybName = hybName, labelVariable = labelVariable, 
    requireLog2=requireLog2)
  qcNames<-as.character(fData(esQC)[, c(reporterGroupName)])
  
  nProbes<-length(probes) 
  for(k in 1:nProbes)
  {
    cat("\n********** k=", k, " *******\n")     
    cat("QC probe=", probes[k], "\n")
    probe.k=match.arg(probes[k], choices=c("biotin", 
      "cy3_hyb", "housekeeping",
      "low_stringency_hyb", "signal", "p95p05"))

    if(probe.k == "signal")
    {
      probs<-c(0.05, 0.25, 0.5, 0.75, 0.95)
      qMat<-apply(datQC, 2, quantile, na.rm=TRUE, prob=probs)

      #if(plotOutPutFlag)
      #{
         outFileName<-paste(projectName, "_", probe.k, "_traj_plot.pdf", sep="")
      #  postscript(outFileName, horizontal=TRUE, paper="letter")
      #} 
      curveNames<-paste("Signal P", round(probs*100,2),sep="")
      plotCurves(qMat, curveNames=curveNames, 
        fileName=outFileName,
        plotOutPutFlag=plotOutPutFlag,
        requireLog2=FALSE, cex=cex,
        ylim=ylim, xlab=xlab, ylab=ylab,
        main=probe.k, lwd=lwd, mar=mar, cex.axis=cex.axis, ...)
      #if(plotOutPutFlag)
      #{
      #  dev.off()
      #}
    } else if (probe.k == "p95p05") {
      qMat<-apply(datQC, 2, quantile, na.rm=TRUE, prob=c(0.05, 0.95))
      ratio<-qMat[2,]/qMat[1,]
      p<-length(ratio)
      if(plotOutPutFlag)
      {
        outFileName<-paste(projectName, "_", probe.k, "_traj_plot.pdf", sep="")
        #postscript(outFileName, horizontal=TRUE, paper="letter")
        pdf(outFileName, width=756, height=576, paper="letter")
      } 
      if(requireLog2)
      { 
        ylab2="log2(p95/p05)"
      } else {
        ylab2="p95/p05"
      }
      par(mar=mar)
      plot(x=1:p, y=ratio, type="l",
        xlab=xlab, ylab=ylab2,
        main="p95/p05", cex=cex,
        ylim=ylim, col=1, lty=1, lwd=lwd, axes=FALSE,
        cex.axis=cex.axis,...)
      box()
      axis(side=2, cex.axis=1, ...)
      axis(side=1, at=1:p, labels=colnames(datQC), las=las, 
        cex.axis=cex.axis, ...)
      par(mar=c(5, 4, 4, 2) + 0.1)
      
      if(plotOutPutFlag)
      {
        dev.off()
      }
    } else {
      pos<-which(qcNames==probe.k)
      len.pos<-length(pos)
      if(len.pos)
      {
        esQC2<-esQC[pos,]
        curveNames<-paste(probe.k, 1:len.pos, sep=" ")
        #if(plotOutPutFlag)
        #{
          outFileName<-paste(projectName, "_", probe.k, "_traj_plot.pdf", sep="")
        #  postscript(outFileName, horizontal=TRUE, paper="letter")
        #} 
        
        datQC2<-getDatColnames(esQC2, hybName = hybName, labelVariable = labelVariable, 
          requireLog2=requireLog2)

        plotCurves(dat=datQC2, curveNames=curveNames,
          fileName = outFileName,
          plotOutPutFlag = plotOutPutFlag,
          requireLog2=FALSE, cex=cex,
          ylim=ylim, xlab=xlab, ylab=ylab, 
          main=probe.k, lwd=lwd, mar=mar, cex.axis=cex.axis, ...)
        #if(plotOutPutFlag)
        #{
        #  dev.off()
        #}
      } else {
        cat("Error! no QC probe", probe.k, " in QC data!\n")
        stop("QC probe not in QC data!\n")
      }
    }
  }
}


# plot trajectories of ratio of 95th percentile to 5th percentile
# of sample probe expression levels across arrays
plotSamplep95p05<-function(es,
    labelVariable="subjID",
    hybName = "Hybridization_Name",
    requireLog2=FALSE,
    projectName="test",
    plotOutPutFlag=FALSE,
    cex=1,
    ylim=NULL,
    xlab="",
    ylab="",
    lwd=1.5,
    mar=c(10, 4, 4, 2) + 0.1,
    las=2,
    cex.axis=1.5,
    title="Trajectory of p95/p05",
    cex.legend=1.5,
    cex.lab=1.5,
    legendPosition="topright",
    cut1=10,
    cut2=6,
    sortFlag = TRUE,
    varSort=c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
    timeFormat=c("%m/%d/%Y", NA, NA),
    verbose=FALSE,
    ...)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }

  if(sortFlag)
  { es<-sortExpressionSet(es, varSort, timeFormat) }

  dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable,
    requireLog2=requireLog2)

  arrayNames.u<-unique(colnames(dat))

  labels<-as.numeric(as.factor(colnames(dat)))
  labels.u<-unique(labels)
  nClusters<-length(labels.u)


  qMat<-apply(dat, 2, quantile, na.rm=TRUE, prob=c(0.05, 0.95))
  ratio<-qMat[2,]/qMat[1,]
  if(verbose)
  {
    cat("sum(ratio<10)>>\n"); print(sum(ratio<10)); cat("\n");
    cat("which(ratio<10)>>\n"); print(which(ratio<10)); cat("\n");
  }

  tmpos<-which(ratio<cut1)
  if(length(tmpos))
  {
    labels2<-labels[tmpos]
    if(verbose)
    { cat("labels for arrays with ratio <", cut1, ">>\n")
      print(table(labels2, useNA="ifany"))
      cat("\n")

      cat("p95 for these arrays>>\n")
      print(qMat[2,tmpos])
      cat("\n")
  
      cat("p05 for these arrays>>\n")
      print(qMat[1,tmpos])
      cat("\n")
  
  
      cat("p95/p05 for these arrays>>\n")
      print(ratio[tmpos])
      cat("\n")
    }
  }

  tmpos<-which(ratio<cut2)
  if(length(tmpos))
  {
    labels2<-labels[tmpos]
    if(verbose)
    {
      cat("labels for arrays with ratio <", cut2, ">>\n")
      print(table(labels2, useNA="ifany"))
      cat("\n")
  
      cat("p95 for these arrays>>\n")
      print(qMat[2,tmpos])
      cat("\n")
  
      cat("p05 for these arrays>>\n")
      print(qMat[1,tmpos])
      cat("\n")
  
  
      cat("p95/p05 for these arrays>>\n")
      print(ratio[tmpos])
      cat("\n")
    }  
  }


  p<-length(ratio)
  if(requireLog2 && is.null(ylab))
  {
    ylab2="log2(p95/p05)"
  } else {
    ylab2="p95/p05"
  }
  if(plotOutPutFlag)
  {
    outFileName<-paste(projectName, "_Samplep95p05_traj_plot.pdf", sep="")
    pdf(outFileName, width=756, height=576, paper="letter")
  }
  par(mar=mar)
  mycols <- grDevices::rainbow(nClusters)
  plot(x=1:p, y=ratio, type="l",
    xlab=xlab, ylab=ylab2,
    main=title, cex=cex,
    ylim=ylim, col=1, pch=labels, lty=1, lwd=lwd, axes=FALSE,
    ...)
  if(labelVariable!="subjID")
  {
    points(x=1:p, y=ratio,
      col=mycols[labels], pch=labels, lwd=lwd,
      ...)
  }
  box()
  axis(side=2, cex.axis=1, ...)
  axis(side=1, at=1:p, labels=colnames(dat), las=las, cex.axis=cex.axis, ...)
  abline(h=cut1, col=1)
  abline(h=cut2, col=1)
  if(labelVariable!="subjID")
  {
    legend(x=legendPosition, legend=arrayNames.u,
      col=unique(mycols[labels]), pch=labels.u,
      cex=cex.legend, bty="o")
  }

  par(mar=c(5, 4, 4, 2) + 0.1)
  if(plotOutPutFlag)
  {
    dev.off()
  }
  res<-list(qMat=qMat, ratio=ratio)
  invisible(res)
}



# plot trajectories of probe expression levels across arrays
plotCurves<-function(dat, curveNames, fileName, 
    plotOutPutFlag=FALSE,
    requireLog2=FALSE, cex=1,
    ylim=NULL, xlab="", ylab="intensity", lwd=3,
    main="Trajectory plot", mar=c(10, 4, 4, 2) + 0.1, las=2, 
    cex.axis=1, ...)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }

  if(requireLog2)
  {
    dat<-log2(dat)
  } 
  if(is.null(ylim))
  { 
    ylim<-range(c(dat), na.rm=TRUE) 
    len.interval<-ylim[2]-ylim[1] 
    ylim[2]<-ylim[2]+len.interval/5
  }
  arrayNames<-colnames(dat)

  n<-nrow(dat)
  p<-ncol(dat)

  if(plotOutPutFlag)
  {
    #postscript(fileName, horizontal=TRUE, paper="letter")
    pdf(fileName, width=756, height=576, paper="letter")
  } 
  mycols <- grDevices::rainbow(n)
  par(mar=mar)
  plot(x=1:p, y=dat[1,], type="l",
    xlab=xlab, ylab=ylab,
    main=main, cex=cex,
    ylim=ylim, col=mycols[1], lty=1, lwd=lwd, axes=FALSE, las=las,
    ...)
  box()
  axis(side=2, cex.axis=1, ...)
  axis(side=1, at=1:p, labels=arrayNames, las=2, cex.axis=cex.axis, ...)
  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x=1:p, y=dat[i,], col=mycols[i], lty=i, lwd=lwd)
    }
  }
  legend("topright", legend=curveNames,
    col=mycols, lty=1:n, cex=cex)
  par(mar=c(5, 4, 4, 2) + 0.1)
  if(plotOutPutFlag)
  {
    dev.off()
  } 
 
}

# plot trajectories of probe expression levels across arrays
quantilePlot<-function(dat, fileName, 
    probs=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
    plotOutPutFlag=FALSE,
    requireLog2=FALSE, 
    sortFlag=TRUE,
    cex=1,
    ylim=NULL, xlab="", ylab="intensity", lwd=3,
    main="Trajectory plot of quantiles", mar=c(15, 4, 4, 2) + 0.1, las=2, 
    cex.axis=1)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }

  if(requireLog2)
  {
    dat<-log2(dat)
  } 
  if(sortFlag)
  {
    # sort columns by MAD
    cat("***** Arrays were sorted by MAD (median absolute deviation)!\n")
    tt.mad<-apply(dat, 2, mad, na.rm=TRUE)
    ttmat<-cbind(tt.mad, 1:length(tt.mad))
    ttmat2<-ttmat[order(ttmat[,1]),]
    dat<-dat[, ttmat2[,2]]
  }
  qMat<-apply(dat, 2, quantile, probs=probs, na.rm=TRUE)
  if(is.null(ylim))
  { 
    ylim<-range(c(qMat), na.rm=TRUE) 
    #len.interval<-ylim[2]-ylim[1] 
    #ylim[2]<-ylim[2]+len.interval/5
  }
  arrayNames<-colnames(dat)
  curveNames<-paste(probs*100, "%", sep="")

  n<-nrow(qMat) # number of quantiles
  p<-ncol(qMat) # number of arrays

  if(plotOutPutFlag)
  {
    #postscript(fileName, horizontal=TRUE, paper="letter")
    pdf(fileName, width=756, height=576, paper="letter")
  } 
  par(mar=mar)
  #mycols <- grDevices::rainbow(n)
  mycols <- 1:n
  plot(x=1:p, y=qMat[1,], type="l",
    xlab=xlab, ylab=ylab,
    main=main, cex=cex,
    ylim=ylim, col=mycols[1], lty=1, lwd=lwd, axes=FALSE, las=las,
    )
  box()
  axis(side=2, cex.axis=1)
  axis(side=1, at=1:p, labels=arrayNames, las=2, cex.axis=cex.axis)
  if(n>1)
  {
    for(i in 2:n)
    {
      lines(x=1:p, y=qMat[i,], col=mycols[i], lty=i, lwd=lwd)
    }
  }
  legend("topright", legend=curveNames,
    col=mycols, lty=1:n, cex=cex)
  par(mar=c(5, 4, 4, 2) + 0.1)
  if(plotOutPutFlag)
  {
    dev.off()
  } 
  invisible(qMat) 
}

# es is an ExpressionSet object
# labelVariable -- the name of a phenotype data variable use to
#                label the arrays in the dendrogram.
getPCAFunc<-function(es, 
  labelVariable="subjID",
  hybName = "Hybridization_Name",
  requireLog2=TRUE,
  corFlag = FALSE
)
{
  dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable,
    requireLog2=requireLog2)

  aa<-apply(dat, 2, function(x) sum(is.na(x)==TRUE))
  pos.NA.arrays<-which(aa==nrow(dat))
  len.NA.arrays=length(pos.NA.arrays)
  if(len.NA.arrays)
  {

    if(len.NA.arrays==ncol(dat))
    {
      stop("The expression levels of all arrays are all missing!\n")
    }
    cat("The expression levels of the following arrays are all missing>>\n")
    print(sampleNames(es)[pos.NA.arrays])
    cat("These arrays will be excluded!\n")

    es<-es[,-pos.NA.arrays] 
    dat<-dat[,-pos.NA.arrays]
  }
  
  aa<-apply(dat, 1, function(x) sum(is.na(x)==TRUE))
  pos.NA.probes <-which(aa==ncol(dat))
  len.NA.probes<-length(pos.NA.probes)
  if(len.NA.probes)
  {

    if(len.NA.probes==nrow(dat))
    {
      stop("The expression levels of all probes are all missing!\n")
    }
    cat("The expression levels of the following probes are all missing>>\n")
    print(featureNames(es)[pos.NA.probes])
    cat("These probes will be excluded!\n")

    es<-es[-pos.NA.probes, ] 
    dat<-dat[-pos.NA.probes,]
  }
 
  # remove probes with missing values
  dat2=na.omit(dat)
  pos=match(rownames(dat2), rownames(dat))
  nNA=nrow(dat)-nrow(dat)
  if(nNA)
  {
    cat("\nWarning: ", nNA, " number of probes containing missing valuse and were deleted\n")
  } 
  
  es2=es[pos,]
  nn<-min(dim(dat2))
  pcs<-prcomp(t(dat2), scale. =corFlag)
  
  res<-list(es.s=es2, pcs=pcs, requireLog2=requireLog2)
  invisible(res)
}

pca2DPlot<-function(pcaObj,
  plot.dim = c(1,2),
  labelVariable="subjID",
  hybName = "Hybridization_Name",
  outFileName="test_pca_raw.pdf",
  title="Scatter plot of pcas",
  plotOutPutFlag=FALSE,
  mar=c(5, 4, 4, 2) + 0.1,
  lwd=1.5,
  equalRange=TRUE,
  xlab=NULL,
  ylab=NULL,
  xlim=NULL,
  ylim=NULL,
  cex.legend=1.5,
  cex=1.5,
  cex.lab=1.5,
  cex.axis=1.5,
  legendPosition="topright",
  ...
)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }

  if(length(plot.dim)!=2 || any(plot.dim<1))
  {
    stop("plot.dim must have exactly 2 positive-integer-value elements!")
  }
  es<-pcaObj$es.s
  if(any(plot.dim > min(dim(es))))
  {
    stop("all elements of plot.dim must be smaller than min(dim(es))!")
  }
  pcs<-pcaObj$pcs
  requireLog2<- pcaObj$requireLog2

  pcaMat<-pcs$x
  aa<-pcaMat[, plot.dim, drop=FALSE]
  egvals<-pcs$sdev^2
  R2vec<-egvals/sum(egvals)

  if(is.null(xlim))
  {
    xlim <- range(aa[,1], na.rm=TRUE)
  }
  if(is.null(ylim))
  {
    ylim <- range(aa[,2], na.rm=TRUE)
  }

  if(equalRange)
  {
    t1<-min(xlim[1], ylim[1])
    t2<-max(xlim[2], ylim[2])
    xlim<-c(t1, t2)
    ylim<-c(t1, t2)
  }


  if(is.null(xlab))
  {
    xlab<-paste("pc", plot.dim[1], " (R2=", round(R2vec[plot.dim[1]],2), ")", sep="")
  }
  if(is.null(ylab))
  {
    ylab<-paste("pc", plot.dim[2], " (R2=", round(R2vec[plot.dim[2]],2), ")", sep="")
  }

  if(labelVariable=="subjID")
  {
    # pca plot
    if(plotOutPutFlag)
    {
      #postscript(outFileName, horizontal=TRUE, paper="letter")
      pdf(outFileName, width=756, height=576, paper="letter")
    }

    par(mar=mar)

    if(is.na(title))
    {
      main=labelVariable
    } else {
      main=title
    }
    plot(x=as.numeric(as.character(aa[,1])),
      y=as.numeric(as.character(aa[,2])),
      xlab=xlab, ylab=ylab,type="p",
      xlim=xlim, ylim=ylim, col=1, pch=1, lwd=lwd,
      cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, ...)
    title(main=main)

    par(mar=c(5, 4, 4, 2) + 0.1)

    if(plotOutPutFlag)
    {
      dev.off()
    }

  } else {
    dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable,
      requireLog2=requireLog2)

    arrayNames.u<-unique(colnames(dat))

    labels<-as.numeric(as.factor(colnames(dat)))
    labels.u<-unique(labels)
    nClusters<-length(labels.u)

    # pca plot
    if(plotOutPutFlag)
    {
      #postscript(outFileName, horizontal=TRUE, paper="letter")
      pdf(outFileName, width=756, height=576, paper="letter")
    }

    par(mar=mar)

    if(is.na(title))
    {
      main=labelVariable
    } else {
      main=title
    }
    mycols <- grDevices::rainbow(nClusters)
    plot(
      x=aa[labels==labels.u[1],1],
      y=aa[labels==labels.u[1],2],
      xlab=xlab, ylab=ylab,type="p",
      xlim=xlim, ylim=ylim, col=mycols[1], pch=labels.u[1], lwd=lwd,
      cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, ...)
    title(main=main)
    if(nClusters>1)
    {
      for(i in 2:nClusters)
      {
        points(x=aa[labels==labels.u[i],1], y=aa[labels==labels.u[i],2],
          col=mycols[i], pch=labels.u[i], lwd=lwd, ...)
      }
    }
    legend(x=legendPosition, legend=arrayNames.u,
      col=mycols, pch=labels.u,
      cex=cex.legend, bty="o")

    par(mar=c(5, 4, 4, 2) + 0.1)

    if(plotOutPutFlag)
    {
      dev.off()
    }
  }
  invisible(pcaMat)
}

pca3DPlot<-function(pcaObj,
  plot.dim=c(1,2,3), 
  labelVariable="subjID",
  hybName = "Hybridization_Name",
  outFileName="test_pca_raw.pdf",
  title="Scatter plot of pcas",
  plotOutPutFlag=FALSE,
  mar=c(5, 4, 4, 2) + 0.1,
  lwd=1.5,
  equalRange=TRUE,
  xlab=NULL,
  ylab=NULL,
  zlab=NULL,
  xlim=NULL,
  ylim=NULL,
  zlim=NULL,
  cex.legend=1.5,
  cex=1.5,
  cex.lab=1.5,
  cex.axis=1.5,
  legendPosition="topright",
  ...
)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }

  if(length(plot.dim)!=3 || any(plot.dim<1))
  {
    stop("plot.dim must have exactly 3 positive-integer-value elements!")
  }
  es<-pcaObj$es.s
  if(any(plot.dim > min(dim(es))))
  {
    stop("all elements of plot.dim must be smaller than min(dim(es))!")
  }
  pcs<-pcaObj$pcs
  requireLog2<- pcaObj$requireLog2

  #pcaMat<-pls::scores(pcs)
  pcaMat<-pcs$x
  aa<-pcaMat[, plot.dim, drop=FALSE]
  #R2vec<-pcs@R2[plot.dim]
  aa<-pcaMat[, plot.dim, drop=FALSE]
  egvals<-pcs$sdev^2
  R2vec<-egvals/sum(egvals)
  
  #R2vec<-pcs@R2[plot.dim]

  if(is.null(xlab))
  {
    xlab<-paste("pc", plot.dim[1], " (R2=", round(R2vec[plot.dim[1]],2), ")", sep="")
  }
  if(is.null(ylab))
  {
    ylab<-paste("pc", plot.dim[2], " (R2=", round(R2vec[plot.dim[2]],2), ")", sep="")
  }
  if(is.null(zlab))
  {
    zlab<-paste("pc", plot.dim[3], " (R2=", round(R2vec[plot.dim[3]],2), ")", sep="")
  }

  if(is.null(xlim))
  {
    xlim <- range(aa[,1], na.rm=TRUE)
  }
  if(is.null(ylim))
  {
    ylim <- range(aa[,2], na.rm=TRUE)
  }
  if(is.null(zlim))
  {
    zlim <- range(aa[,3], na.rm=TRUE)
  }


  if(equalRange)
  {
    t1<-min(c(xlim[1], ylim[1], zlim[1]))
    t2<-max(c(xlim[2], ylim[2], zlim[2]))
    xlim<-c(t1, t2)
    ylim<-c(t1, t2)
    zlim<-c(t1, t2)
  }


  if(labelVariable=="subjID")
  {
    # pca plot
    if(plotOutPutFlag)
    {
      #postscript(outFileName, horizontal=TRUE, paper="letter")
      pdf(outFileName, width=756, height=576, paper="letter")
    }

    par(mar=mar)

    if(is.na(title))
    {
      main=labelVariable
    } else {
      main=title
    }
    scatterplot3d(
      x=as.numeric(as.character(aa[,1])),
      y=as.numeric(as.character(aa[,2])),
      z=as.numeric(as.character(aa[,3])),
      xlab=xlab, ylab=ylab,zlab=zlab,
      xlim=xlim, ylim=ylim, zlim=zlim, 
      color=1, pch=1, lwd=lwd,
      cex.symbols=cex, cex.lab=cex.lab, cex.axis=cex.axis, ...)
    title(main=main)

    par(mar=c(5, 4, 4, 2) + 0.1)

    if(plotOutPutFlag)
    {
      dev.off()
    }

  } else {
    dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable,
      requireLog2=requireLog2)

    arrayNames.u<-unique(colnames(dat))

    labels<-as.numeric(as.factor(colnames(dat)))
    mycols<-labels
    labels.u<-unique(labels)
    nClusters<-length(labels.u)

    # pca plot
    if(plotOutPutFlag)
    {
      #postscript(outFileName, horizontal=TRUE, paper="letter")
      pdf(outFileName, width=756, height=576, paper="letter")
    }

    par(mar=mar)

    if(is.na(title))
    {
      main=labelVariable
    } else {
      main=title
    }
    #mycols <- grDevices::rainbow(nClusters)
    #mycols<-as.numeric(as.factor(colnames(dat)))
    
    scatterplot3d(
      x=aa[,1],
      y=aa[,2],
      z=aa[,3],
      xlab=xlab, ylab=ylab,zlab=zlab, type="p",
      xlim=xlim, ylim=ylim, zlim=zlim,
      color=labels, pch=labels, lwd=lwd,
      cex.symbols=cex, cex.lab=cex.lab, cex.axis=cex.axis, ...)
    title(main=main)
    legend(x=legendPosition, legend=arrayNames.u,
      col=labels.u, pch=labels.u,
      cex=cex.legend, bty="o")

    par(mar=c(5, 4, 4, 2) + 0.1)

    if(plotOutPutFlag)
    {
      dev.off()
    }
  }
  invisible(pcaMat)
}

###########################################################

# get GC samples
separateGCnonGC<-function(es, GCid=c("128115", "Hela", "Brain"))
{
  res<-extractIDs3(pDat=pData(es), mySet=GCid)
  pos.GC<-res$pos
  if(length(pos.GC))
  {
    esGC<-es[, pos.GC]
    es<-es[,-pos.GC]
  } else {
    cat("No GC samples!\n")
    esGC<-NULL
  }

  res<-list(esGC=esGC, esNonGC=es)
  invisible(res)
}




upperTriangle <- function(x, diag=FALSE)
  {
    x[upper.tri(x, diag=diag)]
  }

# es is an ExpressionSet objectthe object for sample probe profiles
# labelVariable -- the name of a phenotype data variable use to
#                label the arrays in the boxplots
# 
# arrayType
R2PlotFunc<-function(es, 
  hybName = "Hybridization_Name",
  arrayType=c("all", "replicates", "GC"),
  GCid = c("128115", "Hela", "Brain"),
  probs = seq(0, 1, 0.25),
  col=gplots::greenred(75),
  labelVariable="subjID",
  outFileName="test_R2_raw.pdf",  
  title="Raw Data R^2 Plot",
  requireLog2=FALSE, 
  plotOutPutFlag=FALSE,
  las=2,
  keysize=1.0,
  margins=c(10,10),
  sortFlag = TRUE,
  varSort=c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
  timeFormat=c("%m/%d/%Y", NA, NA),
  ...
)
{
  if (!interactive())
  {
    options(rgl.useNULL = TRUE)
  }


  arrayType=match.arg(arrayType, choices=c("all", "replicates", "GC"))

  dat<-exprs(es)
  pDat<-pData(es)

  R2vec.within.rep<-NULL
  if(arrayType=="all")
  {
    if(sortFlag)
    { es<-sortExpressionSet(es, varSort, timeFormat) }
    dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable, 
      requireLog2=requireLog2)
    R2Mat<-cor(x=dat, use="complete.obs")^2
  } else if(arrayType=="replicates") {
    # some Subject_ID are NA
    subjIDs.all<-sapply(as.character(pDat[, c(hybName)]), function(x) {
        x<-as.character(x)
        tt<-unlist(strsplit(x, split="_"))
        return(tt[1])
      }
    )
    tt<-table(subjIDs.all, useNA="ifany")
    pos<-which(tt>1)
    subj.r<-names(tt)[pos]
    
    len<-length(subj.r)
    
    pos.rep<-NULL
    for(i in 1:len)
    {
      pos.rep.i<-which(subjIDs.all==subj.r[i])
      pos.rep<-c(pos.rep, pos.rep.i)  
    }
 
    n.rep<-length(pos.rep)
    if(n.rep)
    { 
      es<-es[,pos.rep]
      if(sortFlag)
      { es<-sortExpressionSet(es, varSort, timeFormat) }
      pDat<-pData(es)
      dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable, 
        requireLog2=requireLog2)
      R2Mat<-cor(x=dat, use="complete.obs")^2
 
      # column names
      cn<-as.character(pDat[, c(hybName)])
      # get subject IDs
      subjIDs<-sapply(cn, function(x) {
          x<-as.character(x)
          bb<-unlist(strsplit(x, split="_"))
          return(bb[1])
        }
      )
      # number of unique subjects
      subj.u<-unique(subjIDs)
      n.subj<-length(subj.u)
      for(k in 1:n.subj)
      {
        pos.k<-which(subjIDs==subj.u[k])
        R2Mat.k<-R2Mat[pos.k, pos.k]
        R2vec.k<-c(upperTriangle(R2Mat.k, diag=FALSE))
        R2vec.within.rep<-c(R2vec.within.rep, R2vec.k)
      }
    } else {
      stop("No replicates!\n")
    }
  } else { # GC- genetic control
    es<-extractSamples(es=es, GCflag=TRUE, GCIDsymbols=GCid, hybName=hybName, verbose=FALSE) 
    if(sortFlag)
    { es<-sortExpressionSet(es, varSort, timeFormat) }

    dat<-getDatColnames(es, hybName = hybName, labelVariable = labelVariable, 
      requireLog2=requireLog2)
    R2Mat<-cor(x=dat, use="complete.obs")^2
  }
  rm(dat)
  rm(pDat)

  R2vec<-c(upperTriangle(R2Mat, diag=FALSE))
  qVec<-quantile(R2vec, probs=probs, na.rm=TRUE)
  cat("quantile of R^2>>\n"); print(qVec); cat("\n");
  if(!is.null(R2vec.within.rep))
  { qVec.within.rep<-quantile(R2vec.within.rep, probs=probs, na.rm=TRUE)
    cat("quantile of within-replicate R^2>>\n"); print(qVec.within.rep); cat("\n");
  }

  # draw heatmap
  ##################
  if(plotOutPutFlag)
  {
    pdf(outFileName, width=756, height=576, paper="letter")
    #postscript(outFileName, horizontal=TRUE, paper="letter")
  } 

  gplots::heatmap.2(R2Mat, col=col, Rowv=NULL, Colv=NULL, 
    dendrogram ="none", key=TRUE, keysize=keysize, margins=margins,
    symkey=FALSE, density.info="none", trace="none", main=title, ...)
  
  if(plotOutPutFlag)
  {
    dev.off()
  }

  res<-list(R2Mat=R2Mat, R2vec=R2vec, R2vec.within.rep=R2vec.within.rep)
  invisible(res)
}

getDatColnames<-function(es, hybName="Hybridization_Name", labelVariable = "subjID", requireLog2=FALSE)
{
  dat<-exprs(es)
  pDat<-pData(es)
  if(labelVariable == "subjID")
  {
    colnames(dat)<-as.character(pDat[,c(hybName)])
  } else {
    cn<-colnames(pDat)
    pos<-which(cn == labelVariable)
    if(length(pos))
    {
      colnames(dat)<-as.character(pDat[,labelVariable])
    } else {
      stop("labelVariable is not in the phenotype data!\n")
    }
  }

  if(requireLog2)
  {
    dat<-log2(dat)
  } 

  invisible(dat)

}

####################################
# draw scatter plot for top 20 results
# resFrame is the sorted data frame returned by wrapper functions
#  like lmFitWrapper
scatterPlots=function(
  resFrame, 
  es, 
  col.resFrame=c("probeIDs", "stats", "pval", "p.adj"),
  var.pheno="bmi",
  outcomeFlag = FALSE,
  fitLineFlag = TRUE,
  var.probe="TargetID",
  var.gene="Symbol",
  var.chr="Chr",
  nTop=20,
  myylab="expression level",
  datExtrFunc=exprs,
  fileFlag=FALSE,
  fileFormat="ps",
  fileName="scatterPlots.ps")
{

  fDat=fData(es)
  probeAll=as.character(fDat[, c(var.probe)])
  geneAll=as.character(fDat[, c(var.gene)])
  chrAll=as.character(fDat[, c(var.chr)])

  pDat=pData(es)
  #pheno=as.numeric(as.character(pDat[, c(var.pheno)]))
  pheno=pDat[, c(var.pheno)]

  ####
  frame=resFrame
  frame2=frame[1:nTop,]

  cpg.sel=as.character(frame2[, c(col.resFrame[1])])
  pos=match(cpg.sel, probeAll)
  if(sum(is.na(pos)==TRUE))
  {
    stop("some top probes are not in 'es'!")
  }

  gene.sel=geneAll[pos]
  chr.sel=chrAll[pos]

  dat=datExtrFunc(es)
  dat2=dat[pos,,drop=FALSE]
 
  stats=frame2[, c(col.resFrame[2])]
  pval=frame2[, c(col.resFrame[3])]
  if(length(col.resFrame)==4)
  {
    p.adj=frame[, c(col.resFrame[4])]
  }
 
  if(fileFlag)
  {
    if(fileFormat=="ps")
    {
      postscript(file=fileName, horizontal=TRUE, paper="letter")
    } else if(fileFormat=="pdf"){
      pdf(file=fileName)
    } else if(fileFormat=="jpeg"){
      jpeg(filename=fileName)
    } else {
      png(filename=fileName)
    }
  }

  if(length(col.resFrame)==4)
  {
    for(i in 1:nTop)
    {
      yi=dat2[i,]
      if(outcomeFlag) # i.e. 'pheno' is the outcome variable (y-axis)
      {
        myx=yi
        myy=pheno
        xlab=myylab
        ylab=var.pheno
      } else {
        myx=pheno
        myy=yi
        xlab=var.pheno
        ylab=myylab
      }
      plot(x=myx, y=myy, 
        xlab=xlab, ylab=ylab,
        type="p",
        main=paste("scatterplot of ", var.pheno, 
          "\n", cpg.sel[i], "(", gene.sel[i],
          ", chr=", chr.sel[i], ")", sep=""),
  
        sub=paste("stat=", round(stats[i],2),
           ", pval=", sprintf("%.1e", pval[i]),
           ", p.adj=", sprintf("%.1e", p.adj[i]),sep=""))
      if(fitLineFlag)
      {
        res.lm=lm(myy~myx)
        coef=res.lm$coefficients
        abline(a=res.lm$coefficients[1], b=res.lm$coefficients[2], col=2)
      }
    }
  } else {
    for(i in 1:nTop)
    {
      yi=dat2[i,]
      if(outcomeFlag) # i.e. 'pheno' is the outcome variable (y-axis)
      {
        myx=yi
        myy=pheno
        xlab=myylab
        ylab=var.pheno
      } else {
        myx=pheno
        myy=yi
        xlab=var.pheno
        ylab=myylab
      }

      plot(x=myx, y=myy, 
        xlab=xlab, ylab=ylab,
        type="p",
        main=paste("scatterplot of ", var.pheno, 
          "\n", cpg.sel[i], "(", gene.sel[i],
          ", chr=", chr.sel[i], ")", sep=""),
  
        sub=paste("stat=", round(stats[i],2),
           ", pval=", sprintf("%.1e", pval[i]),
           sep=""))
      if(fitLineFlag)
      {
        res.lm=lm(myy~myx)
        coef=res.lm$coefficients
        abline(a=res.lm$coefficients[1], b=res.lm$coefficients[2], col=2)
      }
    }

  }
  if(fileFlag)
  {
    dev.off()
  }

  invisible(0)
}

# draw parallel boxplots for top 20 results
# resFrame is the sorted data frame returned by wrapper functions
#  like lmFitWrapper
boxPlots=function(
  resFrame, 
  es, 
  col.resFrame = c("probeIDs", "stats", "pval", "p.adj"), 
  var.pheno = "sex", 
  var.probe = "TargetID", 
  var.gene = "Symbol", 
  var.chr = "Chr", 
  nTop = 20, 
  myylab = "expression level", 
  datExtrFunc = exprs, 
  fileFlag = FALSE, 
  fileFormat = "ps", 
  fileName = "boxPlots.ps")
{

  fDat=fData(es)
  probeAll=as.character(fDat[, c(var.probe)])
  geneAll=as.character(fDat[, c(var.gene)])
  chrAll=as.character(fDat[, c(var.chr)])

  pDat=pData(es)
  #pheno=as.numeric(as.character(pDat[, c(var.pheno)]))
  pheno=pDat[, c(var.pheno)]

  ####
  frame=resFrame
  frame2=frame[1:nTop,]

  cpg.sel=as.character(frame2[, c(col.resFrame[1])])
  pos=match(cpg.sel, probeAll)
  if(sum(is.na(pos)==TRUE))
  {
    stop("some top probes are not in 'es'!")
  }

  gene.sel=geneAll[pos]
  chr.sel=chrAll[pos]

  # rows are probes and columns are arrays
  dat=datExtrFunc(es)
  dat2=dat[pos,,drop=FALSE]
 
  stats=frame2[, c(col.resFrame[2])]
  pval=frame2[, c(col.resFrame[3])]
  if(length(col.resFrame)==4)
  {
    p.adj=frame[, c(col.resFrame[4])]
  }
 
  if(fileFlag)
  {
    if(fileFormat=="ps")
    {
      postscript(file=fileName, horizontal=TRUE, paper="letter")
    } else if(fileFormat=="pdf"){
      pdf(file=fileName)
    } else if(fileFormat=="jpeg"){
      jpeg(filename=fileName)
    } else {
      png(filename=fileName)
    }
  }

  if(length(col.resFrame)==4)
  {
    for(i in 1:nTop)
    {
      yi=dat2[i,]
      ttdati=data.frame(yi=yi, pheno=pheno)
      boxplot(yi~pheno,dat=ttdati,
        xlab="", ylab=myylab,
        main=paste("boxplot of ", var.pheno, 
          "\n", cpg.sel[i], "(", gene.sel[i],
          ", chr=", chr.sel[i], ")", sep=""),
        sub=paste("stat=", round(stats[i],2),
           ", pval=", sprintf("%.1e", pval[i]),
           ", p.adj=", sprintf("%.1e", p.adj[i]),sep=""))
  
      stripchart(yi~pheno, dat=ttdati,
              vertical = TRUE, method = "jitter", 
              pch = 21, col = "maroon", bg = "bisque", 
              add = TRUE) 
    }
  } else {
    for(i in 1:nTop)
    {
      yi=dat2[i,]
      ttdati=data.frame(yi=yi, pheno=pheno)
      boxplot(yi~pheno,dat=ttdati,
        xlab="", ylab=myylab,
        main=paste("boxplot of ", var.pheno, 
          "\n", cpg.sel[i], "(", gene.sel[i],
          ", chr=", chr.sel[i], ")", sep=""),
        sub=paste("stat=", round(stats[i],2),
           ", pval=", sprintf("%.1e", pval[i]),
           sep=""))
  
      stripchart(yi~pheno, dat=ttdati,
              vertical = TRUE, method = "jitter", 
              pch = 21, col = "maroon", bg = "bisque", 
              add = TRUE) 
    }

  }
  if(fileFlag)
  {
    dev.off()
  }

  invisible(0)
}

# draw estimated density plots for all arrays
densityPlots=function(
  es, 
  requireLog2 = TRUE,
  myxlab = "expression level", 
  mymain = "density plots",
  datExtrFunc = exprs, 
  fileFlag = FALSE, 
  fileFormat = "ps", 
  fileName = "densityPlots.ps")
{

  nArray=ncol(es)
  dat=datExtrFunc(es)
  if(requireLog2)
  {
    dat=log2(dat)
  }

  dLst=list()
  x=NULL
  y=NULL
  for(i in 1:nArray)
  {
    dLst[[i]]=density(na.omit(dat[,i]))
    x=c(x, dLst[[i]]$x)
    y=c(y, dLst[[i]]$y)
  }
  
  myxlim=range(x, na.rm=TRUE)
  myylim=range(y, na.rm=TRUE)

  if(fileFlag)
  {
    if(fileFormat=="ps")
    {
      postscript(file=fileName, horizontal=TRUE, paper="letter")
    } else if(fileFormat=="pdf"){
      pdf(file=fileName)
    } else if(fileFormat=="jpeg"){
      jpeg(filename=fileName)
    } else {
      png(filename=fileName)
    }
  }

  plot(x=dLst[[1]]$x, y=dLst[[1]]$y, 
    xlim=myxlim, ylim=myylim,
    xlab=myxlab, ylab="Density", type="l",
    main=mymain)
  
  
  for(i in 2:nArray)
  {
    lines(x=dLst[[i]]$x, y=dLst[[i]]$y, 
      xlim=myxlim, ylim=myylim,
      xlab="beta value", ylab="Density")
  }
  
  if(fileFlag)
  {
    dev.off()
  }

  invisible(dLst)
}


