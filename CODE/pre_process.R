
cat(" helooo in pre process")
regulators <- as.matrix(read.table(toString(regulatorGeneNamesFileName )))
targets <- as.matrix(read.table(toString(targetGeneNamesFileName)))
pert_file <- toString(perturbationMatrixFile)
pert_input <- as.matrix((readLines(pert_file))) 



#build perturbation matrix from perturbation list and target gene names when there is no header row
buildPertMatrixWithoutHeader <- function(){

  a = sapply(targets, grepl, pert_input) 
  pert <-t(a)
 
  write.table(pert, "./DATA/YEAST_SUBNETWORK/INPUT/perturbation.txt",row.names=FALSE, col.names=FALSE,quote=FALSE,sep='\t')
  return (pert)
  
}

builPertMatrixWithHeader <-function(headerLine)
{

  
  collist <- firstLineFields
  rowlist <-as.character(unlist(targets))
  a = sapply(headerLine, grepl, pert_input)
  b = sapply(rowlist, grepl, pert_input)
  d=t(b) %*% a
  colnames(d) <- collist
  return (d)
}
#check to see if tdata has header row(condition names)

twoLines <- readLines(targetExpressionFile, n=2)

firstLineFields <- strsplit(twoLines, " ")[[1]]

hdr <- all(is.na(as.numeric(firstLineFields)))
if(!hdr) { 
 
  tdata <- as.matrix(read.table(targetExpressionFile, header=F))
  pert <- buildPertMatrixWithoutHeader()

}else{ 
  tdata <- as.matrix(read.table(targetExpressionFile, header=T))
  pert <- builPertMatrixWithHeader(firstLineFields)

}

#build Regulator Expression File from target Expression File

rdata<- tdata[match(regulators, targets),] 
write.table(rdata, "./DATA/YEAST_SUBNETWORK/INPUt/test/rdata.txt",row.names=FALSE, col.names=FALSE,quote=FALSE,sep='\t')

#building allowed matrix, self-regulation is not allowed

a = sapply(regulators, grepl, targets)
allowed <- t(!a)

