##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions which implement multiple simulations or a single simulation.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# simulatePlatform
# Input: 
#	nSim: Number of simulations
#	seed: Starting seed for simulations
#	platformSimulationParameters: Parameters for setting up the simulations
#	isaSimulationScenarios: Scenarios for response generation in ISAs
#	platformDesign: Parameters defining the platform design, including ISAs
#	entryTime: entryTime for novel ISAs - overwrites the information in platformDesign
#	isaList:	listOfISAs - overwrites the information in platformDesign
#	outPath: 	path were to store simulation data
#	savingGrid:	frequency of storing simulation data in outPath (avoid too large data sets to be moved)
# Output:
#	data.frame with simulation results summarized by function "summarizeSingleSimulation", which is customized function to be defined by user.
simulatePlatform=function(
	nSim,seed,
	platformSimulationParameters,
	isaSimulationScenarios,
	platformDesign,entryTime,isaList,outPath,savingGrid=10){
	platformDesign$nStart=length(entryTime==0)
	platformDesign$entryTime=entryTime
	platformDesign$lateISAList=isaList

	outTMP=list()
	k=1
	k1=1

	scenarioSummary=getScenarioSummary(isaSimulationScenarios,length(entryTime))
	for(i in seq(1,nSim)){
		set.seed(seed+i-1)
		test=runSingleTrialSimulation(platformDesign,
				isaSimulationScenarios,platformSimulationParameters)
		outTMP[[k]]=data.frame(summarizeSingleSimulation(i,test[[1]],test[[2]],test[[3]],scenarioSummary))
		if(k==savingGrid){
			write.table(rbindlist(outTMP),file=paste(outPath,"tmp",k1,seed,"0.csv",sep=""),sep="\t",col.names=TRUE,row.names=FALSE)
			outTMP=list()
			k1=k1+1
			k=0
		}
		k=k+1
	}
	if(nSim>=savingGrid){
		out=list()
		for(i in seq(1,(k1-1))){
		out[[i]]=data.table(read.table(paste(outPath,"tmp",i,seed,"0.csv",sep=""),sep="\t",header=TRUE))
		}
	} else {
		out=outTMP	
	}
	print(paste(outPath,"summaryOut",seed,".csv",sep=""))
	write.table(rbindlist(out),file=paste(outPath,"summaryOut",seed,".csv",sep=""),sep="\t",col.names=TRUE,row.names=FALSE)
	return(rbindlist(out))
}
	
# simulateSingleTrial
# Input: 
#	nSim: Number of simulations
#	seed: Starting seed for simulations
#	platformSimulationParameters: Parameters for setting up the simulations
#	isaSimulationScenarios: Scenarios for response generation in ISAs
#	platformDesign: Parameters defining the platform design, including ISAs
#	entryTime: entryTime for novel ISAs - overwrites the information in platformDesign
#	isaList:	listOfISAs - overwrites the information in platformDesign
# Output:
#	list containing the outputs from a single simulation, in particular:
###		the ISA results and final total data in the first list element
###		the scenario summary in the second list element
###		the allocation rates results passed from the simulation in the third list element.
simulateSingleTrial=function(
	nSim,seed,
	platformSimulationParameters,
	isaSimulationScenarios,
	platformDesign,entryTime,isaList){
	set.seed(seed)
	platformDesign$nStart=length(entryTime==0)
	platformDesign$entryTime=entryTime
	platformDesign$lateISAList=isaList

	outTMP=list()
	k=1
	k1=1
	print(c(1,seed+1-1))
	set.seed(seed+1-1)
	scenarioSummary=getScenarioSummary(isaSimulationScenarios,length(entryTime))
		test=runSingleTrialSimulation(platformDesign,
				isaSimulationScenarios,platformSimulationParameters)
	return(list(test,scenarioSummary,test[[3]]))
	#outTMP=data.frame(summarizeSingleSimulation(1,test[[1]],test[[2]],scenarioSummary))
	#return(outTMP)
}
	
