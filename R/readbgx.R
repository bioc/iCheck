readbgx2 <- function(filename, path=".", sep="\t",
       quote="", header=TRUE, probeStart="[Probes]",
       controlStart="[Controls]", ...) 
{
  tmp = readLines(file.path(path, filename))
  pos = grep("^\\[", tmp)
  tmp.sub = tmp[pos]
  #skip = grep(probeStart, tmp)
  skip = pos[which(tmp.sub==probeStart)]
  #end = grep(controlStart, tmp)
  tmppos = which(tmp.sub==controlStart)
  n.tmp.sub = length(tmp.sub)
  end = pos[tmppos]

  # number of rows of probes
  if(header)
    nrows.p = end-skip-2
  else
    nrows.p = end-skip-1

  if(tmppos<n.tmp.sub)
  { 
    end2 = pos[n.tmp.sub] 
  } else {
    end2 = length(tmp)  
  }
  # number of rows of control probes
  if(header)
    nrows.c =  end2-end-2
  else
    nrows.c =  end2-end-1

  probeAnno = read.delim(file.path(path, filename), sep=sep, header=header, skip=skip, nrows=nrows.p, quote=quote, stringsAsFactors=FALSE, ...)
  controlAnno = read.delim(file.path(path, filename), sep=sep, header=header, skip=end, nrows=nrows.c, stringsAsFactors=FALSE, ...)


#  tmp = probeAnno[1:nrow(controlAnno),]
#  ind = match(colnames(controlAnno), colnames(probeAnno))
#  tmp[,-ind[!is.na(ind)]] = ""
#  tmp[,ind[!is.na(ind)]] = controlAnno[,!is.na(ind)]
#  Status = rep("gene", nrow(probeAnno)+nrow(controlAnno))
#  Status[(nrow(probeAnno)+1):(nrow(probeAnno)+nrow(controlAnno))] = as.character(controlAnno$Reporter_Group_Name)
#  anno = rbind(probeAnno, tmp)
#  anno = cbind(anno,Status)
#  anno

  res<-list(probeAnno=probeAnno, controlAnno=controlAnno)
  invisible(res)
}

readbgx <- function(filename, sep="\t",
       quote="", header=TRUE, probeStart="[Probes]",
       controlStart="[Controls]", ...) 
{
  tmp = readLines(filename)
  pos = grep("^\\[", tmp)
  tmp.sub = tmp[pos]
  #skip = grep(probeStart, tmp)
  skip = pos[which(tmp.sub==probeStart)]
  #end = grep(controlStart, tmp)
  tmppos = which(tmp.sub==controlStart)
  n.tmp.sub = length(tmp.sub)
  end = pos[tmppos]

  # number of rows of probes
  if(header)
    nrows.p = end-skip-2
  else
    nrows.p = end-skip-1

  if(tmppos<n.tmp.sub)
  { 
    end2 = pos[n.tmp.sub] 
  } else {
    end2 = length(tmp)  
  }
  # number of rows of control probes
  if(header)
    nrows.c =  end2-end-2
  else
    nrows.c =  end2-end-1

  probeAnno = read.delim(filename, sep=sep, header=header, skip=skip, nrows=nrows.p, quote=quote, stringsAsFactors=FALSE, ...)
  controlAnno = read.delim(filename, sep=sep, header=header, skip=end, nrows=nrows.c, stringsAsFactors=FALSE, ...)


#  tmp = probeAnno[1:nrow(controlAnno),]
#  ind = match(colnames(controlAnno), colnames(probeAnno))
#  tmp[,-ind[!is.na(ind)]] = ""
#  tmp[,ind[!is.na(ind)]] = controlAnno[,!is.na(ind)]
#  Status = rep("gene", nrow(probeAnno)+nrow(controlAnno))
#  Status[(nrow(probeAnno)+1):(nrow(probeAnno)+nrow(controlAnno))] = as.character(controlAnno$Reporter_Group_Name)
#  anno = rbind(probeAnno, tmp)
#  anno = cbind(anno,Status)
#  anno

  res<-list(probeAnno=probeAnno, controlAnno=controlAnno)
  invisible(res)
}
