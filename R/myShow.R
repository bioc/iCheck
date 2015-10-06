myShow<-function(object, showHistory=FALSE)
{
  if(showHistory)
  { cat('Summary of data information:\n')
    note <- Biobase::notes(object)
    n.names <- names(note)
    for (i in seq(note)) {
            if (!is.null(n.names)) cat('\t', paste(n.names[i], ':\n\t\t', sep=''))
            cat(note[[i]], sep='\n\t\t')
    }
    cat('\nMajor Operation History:\n')
    print(lumi::getHistory(object))
  }
  cat('\nObject Information:\n')
  aa<-as(object, "ExpressionSet")
  print(aa)

  if(showHistory)
  {
    if (!is.null(object@controlData)) {
          cat('Control Data: Available\n')
    } else {
        cat('Control Data: N/A\n')
    }
    cat("QC information: Please run mySummary(x) for details!\n")
  }
}



