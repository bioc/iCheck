###########
# created on June 10, 2010
#  (1) read in data
#
###########

########################################
# subset an ExpressionSet
extractSamples<-function(es,
  GCflag = TRUE,
  GCIDsymbols=c("128115", "Hela", "Brain"), 
  hybName = "Hybridization_Name", 
  verbose=TRUE)
{
  tmp<-extractIDs3(pDat=pData(es), mySet=GCIDsymbols, 
		   hybName=hybName,  
		   verbose=verbose)
  if(is.null(tmp))
  {
    pos<-numeric(0)
  } else {
    pos<-tmp$pos
  }

  if(GCflag)
  {
    if(length(pos))
    {
      es2<-es[, pos]
    } else {
      cat("In function 'extractSamples', no GC arrays in the ExpressionSet! NULL will be returned\n")
      return(NULL)
    }
  } else {
    if(length(pos))
    {
      es2<-es[, -pos]
    } else {
      es2<-es
    }
  }
  return(es2)
}
# obtain Subject ID for all arrays
# and for GC arrays.
# The first part of 'sampleIDs' must be Subject ID
# The first part of 'sampleIDs' and other parts of 'sampleIDs'
#  must be separated by '_'
extractIDs<-function(sampleIDs, 
  mySet=c("128115", "Hela", "Brain"), verbose=TRUE)

{
  hn<-as.character(sampleIDs)
  if(is.null(hn))
  {
    cat("The function argument 'sampleIDs' is empty!\n")
    return(NULL)
  }
  IDs<-sapply(hn, function(x) {
      x<-as.character(x)
      tt<-strsplit(x, split="_")
      res<-tt[[1]][1]
      return(res)
    }
  )
  rm(hn)

  sampleIDs.sub<-numeric(0) 
  pos<-numeric(0)
  if(length(mySet))
  { pos<-which(toupper(IDs) %in% toupper(mySet))
    n.pos<-length(pos)
    if(n.pos)
    {
      sampleIDs.sub<-sampleIDs[pos]
    }   
  } else {
    cat("The function argument 'mySet' is empty!\n")
    return(NULL)
  } 
  res<-list(pos=pos, sampleIDs.sub=sampleIDs.sub, IDs=IDs)
  return(res)
}

extractIDs3<-function(pDat, 
  mySet=c("128115", "Hela", "Brain"), 
  hybName = "Hybridization_Name", 
  verbose = TRUE)

{
  sm<-as.character(pDat[, c(hybName)])

  res<-extractIDs(sampleIDs=sm, mySet=mySet, verbose=FALSE)
  ttpos2<-res$pos

  pos<-unique(ttpos2)

  if(length(pos))
  {
    sampleIDs.sub<-sm[pos]
  } else {
    sampleIDs.sub<-NULL
  }

  hn<-as.character(sm)
  IDs<-sapply(hn, function(x) {
      x<-as.character(x)
      tt<-strsplit(x, split="_")
      res<-tt[[1]][1]
      return(res)
    }
  )
  rm(hn)

  res<-list(pos=pos, sampleIDs.sub=sampleIDs.sub, IDs=IDs)
  return(res)
}



########################################

# chipNames -- a vector of chip names 
# chipDir -- directory of arrays
# metaFile -- file name for meta file
#             meta file should be comma delimited with header line
#             and na.strings="NA"
# metaDir -- directory of meta file
# pDat -- a pheno data frame (rows are subjects and columns are variables)
#         (must contain columns: Hybridization_Name)
# mFileName -- file name including full directory path of manifest file
#         containing manifest data (rows are probes and columns are subjects) 
#          (must contain columns: Array_Address_Id
#          Symbol, Chromosome, and ILMN_Gene)
# GConly -- logical. Indicates if only Genetic Control Samples will be output.
#
# outObjFilePrefix - output object file name without file extension
#
#manifestDir="/proj/reilms/reilm00/GX/Manifests/",
readData<-function(
  path, 
  pDat, 
  pattern="SampleProbeProfile.txt",
  QCpattern="QCProbeProfile.txt",
  GCid = c("128115", "Hela", "Brain"),
  GConly = FALSE, 
  manifestDir="Manifests/",
  nuIDFlag = FALSE,
  combineRule = c("union", "intersection"),
  hybName = "Hybridization_Name",
  chipManifest = "Chip_Manifest",
  probeID = "ProbeID",
  arrayAddressID = "Array_Address_Id", 
  probeSeq = "Probe_Sequence",
  reporterGroupName = "Reporter_Group_Name",
  saveFlag = TRUE, 
  outObjFileDir=NULL,
  outObjFilePrefix="es",
  outObjFilePrefixQC="esQC",
  verbose=TRUE
  )
{
  mywd = getwd()

  # read sample probe profile
  cat("\n*** read sample probe profile>>>\n")
  es.raw<-readData.default(
    path=path, 
    pDat=pDat, 
    pattern=pattern,
    sampleProbeFlag = TRUE,
    GCid = GCid, 
    GConly = GConly, 
    manifestDir=manifestDir,
    nuIDFlag = nuIDFlag,
    combineRule = combineRule,
    hybName = hybName,
    chipManifest = chipManifest,
    probeID = probeID,
    arrayAddressID = arrayAddressID, 
    probeSeq = probeSeq,
    saveFlag = FALSE,
    outObjFileDir=outObjFileDir,
    outObjFilePrefix=outObjFilePrefix,
    verbose=verbose
    )

  # read QC probe profile
  cat("\n*** read QC probe profile>>>\n")
  esQC.raw<-readData.default(
    path=path, 
    pDat=pDat, 
    pattern=QCpattern,
    sampleProbeFlag = FALSE,
    GCid = GCid, 
    GConly = GConly, 
    manifestDir=manifestDir,
    nuIDFlag = nuIDFlag,
    combineRule = combineRule,
    hybName = hybName,
    chipManifest = chipManifest,
    probeID = probeID,
    arrayAddressID = arrayAddressID, 
    probeSeq = probeSeq,
    saveFlag = FALSE,
    outObjFileDir=outObjFileDir,
    outObjFilePrefix=outObjFilePrefixQC,
    verbose=verbose
    )

  sn<-sampleNames(es.raw)
  sn.QC<-sampleNames(esQC.raw)

  pos<-match(sn, sn.QC)
  #if(length(sn) != length(sn.QC) || any(is.na(pos)==TRUE))
  if(!identical(sort(sn), sort(sn.QC)))
  {
    sn.diff<-setdiff(sn, sn.QC)
    snQC.diff<-setdiff(sn.QC, sn)
    if(length(sn.diff))
    {
      cat("\n*************************\n")
      cat("The following arrays which are in sample probe profile file\n")
      cat("but are not in QC probe profile file!\n")
      print(sn.diff) 
      cat("\n")
    }
    if(length(snQC.diff))
    {
      cat("\n*************************\n")
      cat("The following arrays which are in QC probe profile file\n")
      cat("but are not in sample probe profile file!\n")
      print(snQC.diff) 
      cat("\n")
    }

    setwd(mywd)
    stop("array names in sample probe files do not match those in QC probe files")
  } 
  esQC.raw<-esQC.raw[,pos]

  ########
  # added controlData slot to es.raw, and esQC.raw.
  # without controlData slot, we can not do bgAdjust.
  ########
  cat("\n**** adding controlData slot>>\n")
  dat.QC<-exprs(esQC.raw)
  rn<-rownames(dat.QC)
  mat.QC<-as.matrix(dat.QC)

  fDat.QC<-fData(esQC.raw)
  if(!identical(rn, as.character(fDat.QC[, c(probeSeq)])))
  {
    setwd(mywd)
    stop(paste("rownames(dat.QC) != fDat.QC$", probeSeq, sep=""))
  }
  rownames(mat.QC)<-as.character(fDat.QC[, c(reporterGroupName)])
  control<-as.data.frame(mat.QC)
  
  es.raw@controlData<-control
  esQC.raw@controlData<-control

  if(saveFlag)
  {
    if(verbose)
    { cat("\n****************************\n") 
      cat("\nSaving raw data>>>>>>>\n")
    }
  
    if(!is.null(outObjFileDir))
    { 
      setwd(outObjFileDir)
    } else {
      outObjFileDir = getwd()
    }
    outFile.u<-paste(outObjFilePrefix, ".raw.rda", sep="")
    outFileQC.u<-paste(outObjFilePrefixQC, ".raw.rda", sep="")
    if(verbose)
    { cat("output file dir is>>", outObjFileDir, "\n")
      cat("output file is>>", outFile.u, "\n")
      cat("output QC file is>>", outFileQC.u, "\n")
    }
    
    save(es.raw, file=outFile.u)
    save(esQC.raw, file=outFileQC.u)
  
  }
  setwd(mywd)

  res<-list(es.raw=es.raw, esQC.raw=esQC.raw)
  invisible(res)
}

# pDat must contain a column 'Hybridization_Name'
# to uniquely define each array
readData.default<-function(
  path, 
  pDat, 
  pattern="SampleProbeProfile.txt",
  sampleProbeFlag = TRUE,
  GCid = c("128115", "Hela", "Brain"),
  GConly = FALSE, 
  manifestDir="Manifests/",
  nuIDFlag = FALSE,
  combineRule = c("union", "intersection"),
  hybName = "Hybridization_Name",
  chipManifest = "Chip_Manifest",
  probeID = "ProbeID",
  arrayAddressID = "Array_Address_Id", 
  probeSeq = "Probe_Sequence",
  saveFlag = TRUE,
  outObjFileDir=NULL,
  outObjFilePrefix="esGC",
  verbose=TRUE
  )
{
  combineRule = match.arg(arg=combineRule, 
    choices=c("union", "intersection"))

  if(is.null(pDat))
  {
    stop("pDat is NULL!\n")
  }

  # check if the directory specified by
  # 'path' exists
  mywd = getwd()
  setwd(path)
  setwd(mywd)

  # check if the directory specified by
  # 'manifestDir' exists
  setwd(manifestDir)
  setwd(mywd)

  # check if the directory specified by
  # 'outObjFileDir' exists
  if(!is.null(outObjFileDir))
  { 
    setwd(outObjFileDir)
    setwd(mywd)
  } else {
    outObjFileDir <- getwd()
  }

  flag<-requireNamespace("lumi")
  if(!flag)
  { 
    setwd(mywd)
    stop("library 'lumi' is required to do vsn normalization!\n") 
  }
  
  hn<-as.character(pDat[, c(hybName)])
  if(is.null(hn))
  {
    setwd(mywd)
    stop(paste("phenotype data 'pDat' should contain the column titled ", hybName, sep=""))
  }
  if(length(unique(hn)) != length(hn))
  {
    setwd(mywd)
    stop(paste("The value of ", hybName, " must be unique!", sep=""))
  }
  tmpres<-extractIDs3(pDat=pDat, mySet=GCid, 
		      hybName=hybName,  
		      verbose=verbose)
  pos<-tmpres$pos
  n.pos<-length(pos)
  if(verbose)
  {
    cat("\n\n************************\n")
    cat("No. of arrays in pDat=", nrow(pDat), "\n")
    cat("No. of GC arrays in pDat=", n.pos, "\n")
    cat("\n\n************************\n")
  }
  if(GConly)
  {
    if(n.pos)
    {
      pDat<-pDat[pos, ,drop=FALSE] 
    } else {
      setwd(mywd)
      stop("No Genetic Controls samples in the pDat data set!\n")
    }
  }

  setwd(path)

  targs = list.files(path, pattern=pattern, 
    recursive=TRUE, full.names=TRUE)

  nFiles<-length(targs)
  if(verbose)
  {
    cat("\n************\n")
    cat("Number of", pattern, "  files = ", nFiles, "\n")
    cat("targs>>>\n"); print(targs); cat("\n");
    cat("\n************\n")
    if(nFiles==0)
    {
      setwd(mywd)
      stop("No files with naming containing the phrase ", pattern, " are found!\n")
    }

    cat("\n**************************\n")
    cat("Reading ", pattern, " files>>>>\n")
  }
  count <- 0
  for (i in 1:nFiles) 
  {
    if(verbose)
    {
      cat("\n ******* reading files .....\n")
      cat("\ni=", i, " targs[i]>>", targs[i], "\n")
    }
    ########
    # read SampleProbeProfile
    ########
    tmp = lumiR(targs[i], QC=FALSE)

    sn<-sampleNames(tmp)

    # since lumiR will automatically added suffix to duplicated arrays
    # so the following condition should be always be FALSE.
    # If not, then there must be something wrong. We stop the program then
    # to check what's wrong.
    tt<-table(sn)
    posi<-which(tt>1)
    n.posi<-length(posi)
    if(n.posi)
    {
      sn.sub<-names(tt[posi])
      tmpframe<-data.frame(sampleName=sn.sub, fileName=rep(targs[i], n.posi))
      cat("\n************************\n")
      cat("length(unique(sampleNames(tmp))) =", length(unique(sn)) , "< length((sampleNames(tmp))) =", length((sn)) , "\n")
      cat("The following sample ids have duplicates\n")
      print(tmpframe)
      cat("\n************************\n")
      setwd(mywd) 
      stop("sampleNames in a LumiBatch object should be unique!")

    }

    noGC = FALSE

    if(GConly)
    {
      tmpres<-extractIDs(sampleIDs=sampleNames(tmp), mySet=GCid, verbose=verbose)
      pos<-tmpres$pos
      if(length(pos)==0)
      { 
        noGC <- TRUE
      } else {
        tmp<-tmp[,pos]
      } 
    }

    if(!noGC || (!GConly)) # contains GC samples
    {
      count <- count + 1
      if(verbose)
      { cat("dim(tmp)>>\n"); print(dim(tmp)); cat("\n"); 
        if(count > 1)
        {
          cat("\nCombining LumiBatch objects>>>\n")
        }
      }

      # get probe ids
      probeName.i<-as.numeric(as.character(featureNames(tmp)))
  
      # some arrays have different orders of rows
      # so we have to sorted rows for each array
      # sort probe IDs
      tmpmat.i<-cbind(probeName.i, 1:length(probeName.i))
      tmpmat.i2<-tmpmat.i[order(tmpmat.i[,1]),]
      rm(tmpmat.i)
     
      # sort the rows of mat.i by probe IDs
      # the reason to do this is to avoid
      # the cases in which some chips have different orders
      # to put probe IDs (i.e. row names are different by orders)
      x.lumi.i<-tmp[tmpmat.i2[,2],]
      rm(tmp)
      rm(tmpmat.i2)

      # add feature data  
      # match the samples in expression file with those in pData file
      # so that we can know which chip version they from
      # Hence, we can read in feature data file and extract 
      # sequence info for merging purpose
      sampleIDs.i<-sampleNames(x.lumi.i)
      aa<-match(sampleIDs.i, as.character(pDat[, c(hybName)]))
      tmppos<-which(is.na(aa)==FALSE)
      if(length(tmppos))
      { 
        chipMani<-as.character(pDat[, c(chipManifest)])
        chipver<-as.character(chipMani[aa[tmppos]])
        if(any(is.na(chipver)==TRUE))
        {
          setwd(mywd)
          #stop("There are empty chip version in pDat!")
          cat("\nWarning: There are empty chip version in pDat!")
        }
        if(length(table(chipver))!=1)
        { setwd(mywd) 
          print(table(chipver, useNA="ifany"))
          #stop("The arrays in the same SampleProbeProfile are from different chip version!\n")
          cat("\nWarning: The arrays in the same SampleProbeProfile are from different chip version!\n")
        }
        setwd(manifestDir)
        manifestFileVec = list.files(manifestDir, pattern=unique(chipver), 
          recursive=TRUE, full.names=TRUE)
        setwd(mywd)
        if(length(manifestFileVec)!=1)
        { cat("chip version>>\n"); print(chipver); cat("\n"); 
          setwd(mywd)
          #stop("There is no manifest file match the chip version!\n")
          cat("\nWarning: There is no manifest file match the chip version!\n")
        }
      
        manifestFile.i<-manifestFileVec[1] 
        if(count ==1 )
        { manifestFile = manifestFile.i }
        if(count==1 || !identical(manifestFile, manifestFile.i))
        {
          # read feature data file 
          tmp.bgx<-readbgx(filename=manifestFile.i, sep="\t",
            quote="", header=TRUE, probeStart="[Probes]",
            controlStart="[Controls]") 
          if(sampleProbeFlag)
          { 
            probeAnno.i<-tmp.bgx$probeAnno
          } else {
            probeAnno.i<-tmp.bgx$controlAnno
          }
        } 
        if(verbose)
        { cat("Adding manfiest data>>>>\n") }
        tt<-addfDat(x.lumi.i, probeAnno.i, 
          featureID.es=probeID, featureID.fDat=arrayAddressID,
          probeSeq = probeSeq, 
          verbose=verbose)
        if(!is.null(tt))
        {
          x.lumi.i<-tt
        } else {
          setwd(mywd)
          stop("Error occurs in function 'addfDat'!")
        }
      } else {
        setwd(mywd)
        cat("\n********************\n")
        cat("\nSamples in the expression files>>>\n")
        print(sampleIDs.i)
        cat("\n********************\n")
        stop("All samples in expression sets are not in pDat!\n")
      }
      setwd(mywd) 


      # use Probe_Sequence as featureNames
      # but Probe_Sequence might not be unique
      # that is, some Probe_Sequence correspond to multiple ProbeIDs
      # merge duplicates and replace feature names
      #  by featureID.es.
      tt<-mergeDuplicates(x.lumi.i, featureID.es=probeSeq)
      if(!is.null(tt))
      {
        x.lumi.i<-tt
      } else {
        setwd(mywd)
        stop("Error occurred in function 'mergeDuplicates'!")
      }

      manifestFile <- manifestFile.i
      if(count>1)
      { 
        tti<-myCombine(x=x.lumi, y=x.lumi.i, combineRule=combineRule, suffix=paste("_", i, sep=""))
        if(is.null(tti))
        {
          setwd(mywd)
          stop("Error in myCombine!\n")
        }
        x.lumi<-tti
      } else {
        x.lumi <- x.lumi.i
      }
      rm(x.lumi.i)
    } else if(GConly && noGC) {
      if(verbose)
      { cat("\n***********\n")
        setwd(mywd)
        stop("GConly=TRUE, but no GC samples!\n")
        cat("\n***********\n")
      }
    }
  }
  setwd(mywd)
  if(count==0)
  {
    setwd(mywd)
    stop("No required expression data avaiable!\n")
  }

  if(verbose)
  { cat("\n****************************\n") 
    cat("merging expression data with phenotype data>>>\n")
  }

  newname.mRNA=as.character(sampleNames(x.lumi))
  subj.pDat<-as.character(pDat[, c(hybName)])

  tt<-match(as.character(newname.mRNA), as.character(subj.pDat))
  pos<-which(is.na(tt)==FALSE)
  if(length(pos)==0)
  {
    setwd(mywd)
    stop("subjects in pDat do not match those in expression data!\n")
  }
  if(length(pos)<length(newname.mRNA))
  {
    if(verbose)
    { cat("\n*********\n")
      cat("\nWarnings: some subjects in expression data do not match those in pDat!\n")
    }
    pos.na<-which(is.na(tt)==TRUE)
    if(verbose)
    {
      cat("These subjects are>>\n")
      print(newname.mRNA[pos.na])
      cat("\n")
      cat("These arrays will be dropped!\n")
      cat("\n*********\n")
    }
  }
  x.lumi<-x.lumi[,pos]
  newname.mRNA<-sampleNames(x.lumi)
  pDat<-pDat[tt[pos],]
  subj.pDat<-as.character(pDat[, c(hybName)])

  flag<-identical(subj.pDat, newname.mRNA)
  if(!flag)
  { 
    setwd(mywd)
    stop("subjIDs in SampleProbe are different from those in phenoData!\n") 
  }

  rownames(pDat)=as.character(pDat[,c(hybName)])
  Biobase::sampleNames(x.lumi)=as.character(pDat[,c(hybName)])

  Biobase::pData(x.lumi)<-pDat
  es<-x.lumi
  rm(x.lumi)
  rm(pDat)

  pDat<-pData(es)

  # chip version
  if(verbose)
  { cat("\n************\n")
    cat("Frequencies of Chip versions>>\n")
    print(table(as.character(pDat[, c(chipManifest)]), useNA="ifany"))
    cat("\n")
    cat("\n************\n")
  }

#  # sort the arrays according
#  #   Batch_Run_Date
#  #   chip barcode (Chip_Barcode)
#  #   array position (Chip_Address)
#  if(verbose)
#  { cat("\n****************************\n") 
#    cat("\n Sort arrays according to Batch_Run_Date, Chip_Barcode, and Chip_Address>>>\n")
#  }
#  es.raw<-sortExpressionSet(es)

  es.raw<-es
  rm(es)

  # get nuIDs
  if(nuIDFlag)
  {
    fDat<-fData(es.raw)
    nuIDs<-seq2id(as.character(fDat[, c(probeSeq)]))
    fDat$nuIDs<-nuIDs
    Biobase::fData(es.raw)<-fDat
    #featureNames(es.raw)<-nuIDs
  }

  if(verbose)
  { cat("dim(es.raw)>>\n"); print(dim(es.raw)); cat("\n"); }

  if(verbose)
  { cat("\n**************************\n")
    cat("Dimension of raw LumiBatch>>>>\n")
    print(dim(es.raw))
    cat("\n")
  }
  if(saveFlag)
  { 

    if(verbose)
    { cat("\n****************************\n") 
      cat("\nSaving raw data>>>>>>>\n")
    }

    if(!is.null(outObjFileDir))
    { 
      setwd(outObjFileDir)
    } else {
      outObjFileDir = getwd()
    }
    outFile.u<-paste(outObjFilePrefix, ".raw.rda", sep="")
    if(verbose)
    { cat("output file dir is>>", outObjFileDir, "\n")
      cat("output file is>>", outFile.u, "\n")
    }
  
    save(es.raw, file=outFile.u)

  }
  setwd(mywd)

  invisible(es.raw)
}

# 
# fixed a bug on Jan 18, 2013
# fDat must have column: Array_Address_Id
# fDat must match the version of chip creating es 
addfDat<-function(es, fDat, featureID.es="ProbeID",
  featureID.fDat="Array_Address_Id", probeSeq = "Probe_Sequence", 
  verbose=TRUE)
{
  if(is.null(fDat[, c(featureID.fDat)]))
  {
    cat("Error: fDat does not have column titled", featureID.fDat, "!\n")
    return(NULL) 
  }
  if(!identical(featureNames(es), fData(es)[, c(featureID.es)]))
  {
    cat("Error: featureNames(es) not equal to ", featureID.es, "!\n") 
    return(NULL)
  }
  if(verbose)
  { cat("\n****************************\n") 
    cat("Adding feature data >>>>\n")
  }
  # update feature data
  ff<-featureNames(es)
  ff.fDat<-as.character(fDat[, c(featureID.fDat)])

  aa<-intersect(ff, ff.fDat) 
  if(length(aa))
  {
    fDat.es<-Biobase::fData(es)
    # check if there are any overlapping column names
    # except featureID.es and featureD.fDat
    int<-intersect(colnames(fDat.es), colnames(fDat))
    len<-length(int)
    if(len)
    {
      int2<-setdiff(int, c(featureID.es, featureID.fDat))
      if(length(int2))
      {
        fDat<-fDat[, -which(colnames(fDat) %in% int2), drop=FALSE]
      }
    }
    #pos.sel<-match(ff, ff.fDat)
    #fDat2<-fDat[pos.sel[which(is.na(pos.sel)==FALSE)],]
    bb2<-merge(x=fDat.es, y=fDat, by.x=featureID.es, 
        by.y=featureID.fDat, all.x=TRUE, sort=FALSE)

    rm(fDat.es)
    #rm(fDat)

    # update the feature data of 'es'
    pos<-match(ff, as.character(bb2[, c(featureID.es)]))
    bb3<-bb2[pos,]
    #ff2<-ff[pos]
    rm(bb2)

    if(!identical(as.character(bb3[, c(featureID.es)]), ff))
    {
      cat("Error: bb3[, ", featureID.es, "] != ff!\n")
      return(NULL)
    }
    
    #rownames(bb3)<-ff2
    rownames(bb3)<-ff
    Biobase::fData(es)<-bb3

    pos.del<-which(is.na(as.character(bb3[, c(probeSeq)]))==TRUE)
    if(length(pos.del))
    {
      cat("Warning: the following probes have missing", probeSeq, " and will be deleted>>>\n")
      print(fData(es)[pos.del,])
      cat("\n")
      es<-es[-pos.del,]
    }
    invisible(es)
  
  } else {
    cat("Error: No overlapping between featureNames(es) and ", featureID.fDat,"!\n")
    return(NULL)
  }
}



# Sort arrays according to Batch_Run_Date, Chip_Barcode, and Chip_Address
sortExpressionSet<-function(es, 
  varSort = c("Batch_Run_Date", "Chip_Barcode", "Chip_Address"), 
  timeFormat = c("%m/%d/%Y", NA, NA)
)
{
  pDat<-pData(es)
  cn<-colnames(pDat)
 
  len<-length(varSort)
  if(len!=length(timeFormat))
  { stop("varSort should have the same length as timeFormt!") }

  pos.nonNA<-which(is.na(timeFormat)==FALSE)
  len.t<-length(pos.nonNA)

  if(len.t)
  {
    for(i in 1:len.t)
    {
      ttpos<-which(cn==varSort[pos.nonNA[i]])
      tt<-attributes(pDat[, ttpos])$class[1]
      tt2<-substr(tt, start=1, stop=4)
      tt3<-as.character(pDat[, ttpos])
      if(!is.null(tt) && tt2=="POSI")
      {
  
      } else {
        # %Y not %y
        # %Y Year with century
        # %y Year without century (00-99)
        tt4<-strptime(tt3, format=timeFormat[pos.nonNA[i]])
        pDat[, ttpos]<-as.character(tt4)
      }
    }  
  }  

  pDat$pos<-1:nrow(pDat)

  if(len>1)
  { 
    pDat2<-pDat[do.call(order, pDat[, c(varSort)]),]
  } else {
    pDat2<-pDat[order(pDat[, c(varSort)]),]
  }

  pos.s<-as.numeric(as.character(pDat2$pos))
  es<-es[,pos.s]
  rm(pDat2)
  invisible(es)
}


# merge duplicates and replace feature names
#  by featureID.es.
# copy and revised from function lumi::lumiR
# and from method Biobase::setReplaceMethod
# (1) To simplify the process, we use sample average for duplicated probes
# for 'exprs', 'se.exprs', 'detection', etc.
# 'lumiR' uses more complicated formula to merge duplicated probes
#  for 'se.exprs', 'detection', etc.  
# (2) 'lumiR' uses 'probeID' as id. But some time, a Probe_Sequence 
#     can have 2 different probeIDs
mergeDuplicates<-function(x.lumi.i, featureID.es="Probe_Sequence")
{
  try(id<-as.character(fData(x.lumi.i)[, featureID.es]), silent=TRUE)
  if(length(id)==0)
  {
    cat("Error:", featureID.es, " is not in fData(x.lumi.i)!\n")
    return(NULL)
  }
  dupId <- unique(id[duplicated(id)])

  obj<-x.lumi.i
  if (length(dupId) > 0) {
      cat("Duplicated Probe_Sequence found and were merged!\n")

      storage.mode <- Biobase::storageMode(obj)
      switch(storage.mode,
             "lockedEnvironment" = {
                 aData <- Biobase::copyEnv(obj@assayData)

                 eltVec<-Biobase::assayDataElementNames(obj)
                 for(elt in eltVec)
                 {
                    rmInd <- NULL
                    for (dupId.i in dupId) {
                      selInd.i <- which(id == dupId.i)
                      aData[[elt]][selInd.i[1], ] <- colMeans(aData[[elt]][selInd.i, , drop = FALSE], na.rm=TRUE)
                      rmInd <- c(rmInd, selInd.i[-1])
                    }
                    if(length(rmInd))
                    { 
                      aData[[elt]]<-aData[[elt]][-rmInd,, drop=FALSE] 
                    }
                 }
               #  Biobase:::assayDataEnvLock(aData)
                 obj@assayData <- aData
             },
             "environment" = {
                 aData <- obj@assayData
                 eltVec<-Biobase::assayDataElementNames(obj)
                 for(elt in eltVec)
                 {
                    rmInd <- NULL
                    for (dupId.i in dupId) {
                      selInd.i <- which(id == dupId.i)
                      aData[[elt]][selInd.i[1], ] <- colMeans(aData[[elt]][selInd.i, , drop = FALSE], na.rm=TRUE)
                      rmInd <- c(rmInd, selInd.i[-1])
                    }
                    if(length(rmInd))
                    {
                      aData[[elt]]<-aData[[elt]][-rmInd,, drop=FALSE]
                    }
                 }
                 obj@assayData <- aData
             }
      )
      #obj
      if (!is.null(fData(obj)) && length(rmInd)) 
      {
        Biobase::fData(obj)<-fData(obj)[-rmInd,,drop=FALSE]
      }
  }
  Biobase::featureNames(obj)<-as.character(fData(obj)[, c(featureID.es)])

  invisible(obj)
}


