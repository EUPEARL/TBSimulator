##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file sources the generic platform simulator source code.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

library(randomizeR)
library(data.table)
library(mvtnorm)

source(paste(pathLib,"allocationFunctions.r",sep=""))
source(paste(pathLib,"dataGenerationFunctions.R",sep=""))
source(paste(pathLib,"dataProcessingFunctions.R",sep=""))
source(paste(pathLib,"genericAnalysisFunctions.r",sep=""))
source(paste(pathLib,"mainSimulationWrapper.r",sep=""))
source(paste(pathLib,"monitoringFunctions.R",sep=""))
source(paste(pathLib,"simulationFunctions.r",sep=""))
source(paste(pathLib,"supportTools.r",sep=""))