# es is a LumiBatch object
LumiBatch2Table<-function(
  es, 
  probeID.var="ProbeID",
  gene.var="Symbol",
  chr.var="Chromosome",
  sep=",", 
  quote=FALSE,
  filePrefix="test", 
  fileExt="csv")
{
  dat<-exprs(es)
  pDat<-pData(es)
  fDat<-fData(es)
  detect<-lumi::detection(es)
  
  idMat<-fDat[,c(probeID.var, gene.var, chr.var)]
  
  dat<-cbind(idMat, dat)
  colnames(dat)[1]="probeID"
  
  
  fileName.dat<-paste(filePrefix, "_exprs.", fileExt, sep="")
  write.table(dat, file=fileName.dat, sep=sep, row.names=FALSE, col.names=TRUE, quote=quote)
  
  if( !is.null(detect))
  {
    detect<-cbind(idMat, detect)
    colnames(detect)[1]="probeID"
    fileName.detect<-paste(filePrefix, "_detection.", fileExt, sep="")
    write.table(detect, file=fileName.detect, sep=sep, row.names=FALSE, col.names=TRUE, quote=quote)
  }
  
  fileName.pdat<-paste(filePrefix, "_pDat.", fileExt, sep="")
  write.table(pDat, file=fileName.pdat, sep=sep, row.names=FALSE, col.names=TRUE, quote=quote)
  
  # some times fDat=NULL
  if(!is.null(fDat))
  {
    fileName.fdat<-paste(filePrefix, "_fDat.", fileExt, sep="")
    write.table(fDat, file=fileName.fdat, sep=sep, row.names=FALSE, col.names=TRUE, quote=quote)
  }
}