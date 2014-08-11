#args <- commandArgs(trailingOnly = TRUE)
cat("hellooo in run_netprophet")
args <- c(
  "./DATA/YEAST_SUBNETWORK/INPUT/data.expr",      
  "./DATA/YEAST_SUBNETWORK/INPUT/rdata.expr", 
  "./DATA/YEAST_SUBNETWORK/INPUT/allowed.adj",
  "./DATA/YEAST_SUBNETWORK/INPUT/data.pert.adj",
  "./DATA/YEAST_SUBNETWORK/INPUT/signed.desig.adj",
  "1",
  "1",
  "lasso.adjmtr",
  "combined_model.adjmtr",
  "./DATA/YEAST_SUBNETWORK/OUTPUT/GLOBAL_SHRINKAGE",
  "combined_model.adjlst",
  "./DATA/YEAST_SUBNETWORK/INPUT/tf.orfs",
  "./DATA/YEAST_SUBNETWORK/INPUT/orfs"
        
)


 
 microarrayFlag <- as.integer(args[6])
 nonGlobalShrinkageFlag <- as.integer(args[7])
 lassoAdjMtrFileName <- toString(args[8])
 combinedAdjMtrFileName <- toString(args[9])
 outputDirectory <- tempdir
 combinedAdjLstFileName <- toString(args[11])
 

#targetExpressionFile <- toString(args[1])
#regulatorExpressionFile <- toString(args[2])
#allowedMatrixFile <- toString(args[3])
#perturbationMatrixFile <- toString(args[4])
# differentialExpressionMatrixFile <- toString(args[5])
# microarrayFlag <- as.integer(args[6])
# nonGlobalShrinkageFlag <- as.integer(args[7])
# lassoAdjMtrFileName <- toString(args[8])
# combinedAdjMtrFileName <- toString(args[9])
# outputDirectory <- toString(args[10])
# combinedAdjLstFileName <- toString(args[11])
# regulatorGeneNamesFileName <- toString(args[12])
# targetGeneNamesFileName <- toString(args[13])


#   
 tdata <- inputTables$tdata
 rdata <- inputTables$rdata
 pert <-inputTables$pert
 allowed <- inputTables$allowed
 de_component <- inputTables$de_component
 regulatorsNames <- inputTables$regulatorsName
  targetsNames <- inputTables$targetsNames
unchanged_de_component <- inputTables$de_component

# rdata <- as.matrix(read.table(regulatorExpressionFile))
# allowed <- as.matrix(read.table(allowedMatrixFile))
# pert <- as.matrix(read.table(perturbationMatrixFile))

targets <- seq(dim(tdata)[1])

source("./CODE/global.lars.regulators.r", local=TRUE)


if(microarrayFlag == 0) {
	##RNA-Seq Data
	tdata <- log(tdata+1)/log(2)
	rdata <- log(rdata+1)/log(2)
}

## Center data
tdata <- tdata - apply(tdata,1,mean)

## Scale data
t.sd <- apply(tdata,1,sd)
t.sdfloor <- mean(t.sd) + sd(t.sd)
t.norm <- apply(rbind(rep(t.sdfloor,times=length(t.sd)),t.sd),2,max) / sqrt(dim(tdata)[2])
tdata <- tdata / ( t.norm * sqrt(dim(tdata)[2]-1) )
#
r.sd <- apply(rdata,1,sd)
r.sdfloor <- mean(r.sd) + sd(r.sd)
r.norm <- apply(rbind(rep(r.sdfloor,times=length(r.sd)),r.sd),2,max) / sqrt(dim(rdata)[2])
rdata <- rdata / ( r.norm * sqrt(dim(rdata)[2]-1) )

## Compute unweighted solution
prior <- matrix(1,ncol=dim(tdata)[1] ,nrow=dim(rdata)[1] )

## TODO: Both seed and # of cv folds (in global.lars.regulators.r) are parameters that should be exposed to the user
seed <- 747
set.seed(seed)

if (nonGlobalShrinkageFlag == 1) {
	uniform.solution <- lars.local(tdata,rdata,pert,prior,allowed)
} else {
	uniform.solution <- lars.multi.optimize(tdata,rdata,pert,prior,allowed)
}

lasso_component <- uniform.solution[[1]]

results$lasso_component <- (lasso=as.table(lasso_component))
#write.table(lasso_component,file.path(outputDirectory,lassoAdjMtrFileName),row.names=FALSE,col.names=FALSE,quote=FALSE)

## Perform model averaging to get final NetProphet Predictions
source("./CODE/combine_models.r", local=TRUE)

# if( file.exists(inputFiles$regulatorGeneNamesFileName) & file.exists(inputFiles$targetGeneNamesFileName)){
# 	source("./CODE/make_adjacency_list.r", local=TRUE)
# }
if( !is.null(regulatorsNames) & !is.null(targetsNames)){
  source("./CODE/make_adjacency_list.r", local=TRUE)
}

