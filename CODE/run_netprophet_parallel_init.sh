#!/bin/bash
targetExpressionFile=${1}
regulatorExpressionFile=${2}
allowedMatrixFile=${3}
perturbationMatrixFile=${4}
differentialExpressionMatrixFile=${5}
microarrayFlag=${6}
nonGlobalShrinkageFlag=${7}
lassoAdjMtrFileName=${8}
combinedModelAdjMtrFileName=${9}
outputDirectory=${10}
combinedAdjLstFileName=${11}
regulatorGeneNamesFileName=${12}
targetGeneNamesFileName=${13}


echo "calling mpirun now"
mpirun -np 11 R --no-save -q --args ${targetExpressionFile} ${regulatorExpressionFile} ${allowedMatrixFile} ${perturbationMatrixFile} ${differentialExpressionMatrixFile} ${microarrayFlag} ${nonGlobalShrinkageFlag} ${lassoAdjMtrFileName} ${combinedModelAdjMtrFileName} ${outputDirectory} ${combinedAdjLstFileName} ${regulatorGeneNamesFileName} ${targetGeneNamesFileName} < run_netprophet_parallel_init.r > ${outputDirectory}/r.out

