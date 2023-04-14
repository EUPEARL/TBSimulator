##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions which execute the main simulation of one single platform trial.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# getRecruitmentRate
# Input: 
#	t: Time in platform
#	recruitmentPars: Vector of length 3 with minimum rate, maximum rate and ED50 of rate
# Output:
#	intensity for simulation of Poisson distributed entry times.
getRecruitmentRate=function(t,recruitmentPars){
	out=recruitmentPars[1]+recruitmentPars[2]*t/(recruitmentPars[3]+t)
}

# recruitment
# Input: 
#	t: Time in platform
#	recruitmentPars: Vector of length 3/4
#	recruitmentType: if not separate, the EMAX type will be used with 3 parameters. If separate, a step function will be used to define the recruitment rate.
# Output:
#	Number of subjects recruited in time t.
recruitment=function(t,recruitmentPars,recruitmentType){
	if(recruitmentType!="separate"){
		out=rpois(1,lambda=getRecruitmentRate(t,recruitmentPars))
	} else {
		if(t<recruitmentPars[2]){
			out=rpois(1,lambda=recruitmentPars[1])	
		} else if(t<recruitmentPars[4]){
			out=rpois(1,lambda=recruitmentPars[1]+recruitmentPars[3])	
		} else {
			out=rpois(1,lambda=recruitmentPars[3])	
		}
	}
	return(out)
}

# getISASimulationParameters
# Input: 
#	design
#	isaSimulationScenarios
# Output:
#	selected effect assumptions for the isa with submitted design.
getISASimulationParameters=function(design,isaSimulationScenarios){
	out=isaSimulationScenarios$scenarioList[[design$effectScenario[[1]]]]
	return(out)
}

# updateSimulationParameters
# Input: 
#	isaDesignList: The list of all currently open ISAs
#	isaSimulationParametersIn: The current list of scenarios for currently open ISAs
#	isaSimulationScenarios: The list from which new scenarios will be selected.
# Output:
#	isaSimulationParameters: List with the scenarios for the designs.
updateSimulationParameters=function(isaDesignList,isaSimulationParametersIn,isaSimulationScenarios,platformDesign){
	if(!is.null(isaSimulationParametersIn)){
		isaSimulationParameters=isaSimulationParametersIn
		nISA=length(isaDesignList)
		nScen=length(isaSimulationParameters)
	} else {
		isaSimulationParameters=list()
		nISA=length(isaDesignList)
		nScen=0
	}

	if(nISA>nScen){
		for(i in seq(nScen+1,nISA)){
			if(i<=2){
				isaSimulationParameters[[i]]=getISASimulationParameters(isaDesignList[[i]],isaSimulationScenarios)
			}
			else {
				if(!isaSimulationScenarios$randomScenario){
					isaSimulationParameters[[i]]=getISASimulationParameters(isaDesignList[[i]],isaSimulationScenarios)
				} else {
					tmp=runif(1)
					index=min(which(cumsum(isaSimulationScenarios$scenarioProb)>tmp))
					isaSimulationParameters[[i]]=isaSimulationScenarios$scenarioList[[index]]
				}
			}
		}
	}		
	return(isaSimulationParameters)
}

# getInitialDesignList
# Input: 
#	platformDesign: The building instructions for the platform
# Output:
#	isaDesignList: Initial design list, with which the platform is starting
getInitialDesignList=function(platformDesign){
	nStart=sum(platformDesign$entryTime==0)
	isaDesignList=list()
	for(i in seq(1,nStart)){
		isaDesignList[[i]]=platformDesign$lateISAList[[i]]
	}
	return(isaDesignList)
}

# updateISAList
# Input: 
#	isaDesignList: The list of ISAs opened so far in the platform
#	dataSummary: A table containing information on sample size for different ISAs and arms within ISAs (if multi-armed)
# Output:
#	isaDesignList: Updated ISA design list, i.e. with new status, current size and control size.
updateISAList=function(isaDesignList,dataSummary){
	nISA=length(isaDesignList)
	controlData=dataSummary[5,2]
	k=1
	newControlSubjects=controlData-isaDesignList[[2]]$currentSize
	for(i in seq(1,nISA)){
		
		nArms=length(isaDesignList[[i]]$targetSize)
		subData=dataSummary[,(k:(k+nArms-1))]
		
		### Here only concurrent subjects are counted.
		if(sum(isaDesignList[[i]]$currentSize)<sum(isaDesignList[[i]]$targetSize)){
			isaDesignList[[i]]$controlSize=isaDesignList[[i]]$controlSize+newControlSubjects
		}

		for(j in seq(1,nArms)){
			isaDesignList[[i]]$currentSize[j]=dataSummary[5,k+(j-1)]
		}
		
		if(nArms==1){
			if(subData[5]>0){
				startDate=subData[3]
				endDate=subData[4]
			} else {
				startDate=0
				endDate=0
			}
		} else if(sum(subData[5,]>0)){
			startDate=min(subData[3,])
			endDate=max(subData[4,])
		} else {
			startDate=0
			endDate=0
		}
	
		if(sum(isaDesignList[[i]]$currentSize)>=sum(isaDesignList[[i]]$targetSize)){
			isaDesignList[[i]]$activeISA=FALSE
		}

		k=k+nArms
	}

	openEnrollmentPerISA=getVectorFromList(isaDesignList,"activeISA",sumElements=TRUE)

	if(any(openEnrollmentPerISA[3:nISA]==1)){
		isaDesignList[[2]]$activeISA=TRUE
	} else {
		isaDesignList[[2]]$activeISA=FALSE
	}
	
	return(isaDesignList)
}

# addNewISA
# Input: 
#	time: time in platform
#	isaDesignListTMP: Current flat isa design list (fewer levels of lists)
#	isaDesignList: Original deep isa design list (many list levels)
#	platformDesign: design instructions, such as entry time for platform.
# Output:
#	isaDesignList: Updated deep isaDesignList
#	new: numberOfAddedISAs
addNewISA=function(time,isaDesignListTMP,isaDesignList,platformDesign){
	nISA=length(isaDesignList)
	new=0
	if(any(platformDesign$entryTime==time)){
		isaDesignListTMP[[2]]$activeISA=TRUE
		isaDesignList=getDeepStructure(isaDesignListTMP,isaDesignList)
		nNewISA=length(which(platformDesign$entryTime==time))
		for(j in seq(1,nNewISA)){
			if(platformDesign$fixedISAOrder){
				index=which(platformDesign$entryTime==time)[j]
			} else {
				x=runif(1)
				index=min(which(cumsum(platformDesign$isaProb)>x))
			}
			isaDesignList[[nISA+j]]=platformDesign$lateISAList[[index]]
		}
		new=nNewISA
	}
	return(list(isaDesignList,new))
}

# getFlatStructure
# Input: 
#	designList: ISAList, which contains lists of lists with parameters for ISAs
#	simulationParameterList: List containing the simulation parameters for response generation in ISA
# Output:
#	outList: An ISAList, which has simulationParameters added and a flat structure
getFlatStructure=function(designList,simulationParameterList){
	nISA=length(designList)
	outList=list()
	if(nISA>0){
		for(i in seq(1,nISA)){
			outList[[i]]=list()
			l=1
			nElementsISA=length(designList[[i]])
			for(j in seq(1,nElementsISA)){
				nElementsDesign=length(designList[[i]][[j]])
				for(k in seq(1,nElementsDesign)){
					outList[[i]][[l]]=designList[[i]][[j]][[k]]
					names(outList[[i]])[l]=names(designList[[i]][[j]])[k]
					l=l+1
				}
			}
			outList[[i]][[l]]=simulationParameterList[[i]]
			names(outList[[i]])[l]="endpointParameters"
		}
	}
	outList[[2]]$activeISA=TRUE
	return(outList)
}		

# getDeepStructure
# Input: 
#	flatDesignList: ISAList which is in flat structure
#	originalDesignList: Original ISAList in deep structure
# Output:
#	outList: An ISAList, which is again in deep structure, but contains elements from flatDesignList
getDeepStructure=function(flatDesignList,originalDesignList){
	nISA=length(flatDesignList)
	outList=originalDesignList
	if(nISA>0){
		for(i in seq(1,nISA)){
			l=1
			nElementsISA=length(originalDesignList[[i]])
			for(j in seq(1,nElementsISA)){
				nElementsDesign=length(originalDesignList[[i]][[j]])
				for(k in seq(1,nElementsDesign)){
					outList[[i]][[j]][[k]]=flatDesignList[[i]][[l]]
					l=l+1
				}
			}
		}
	}
	return(outList)
}		

# prepareISADecisionResult
# Input: 
#	time: time in platform
#	isaDesignList: isaDesignList
# Output:
#	isaDesignList: isaDesignList with added placeholders for decisions, interimPassed,analysisTime, current size and control size
#		The getDecisionMat is a function to be specified in the custom functions.
prepareISADecisionResult=function(time,isaDesignList){
	nISA=length(isaDesignList)
	if(nISA>2){
	for(i in seq(3,nISA)){
		if(is.na(isaDesignList[[i]]$dashboard$decisions[1])){
			nInterim=length(isaDesignList[[i]]$decisionMaking$interimCutoff)
			nArms=length(isaDesignList[[i]]$designElements$targetSize)
			isaDesignList[[i]]$dashboard$decisions=getDecisionMat(isaDesignList[[i]])
			isaDesignList[[i]]$dashboard$interimPassed=rep(0,times=nInterim)
			isaDesignList[[i]]$dashboard$analysisTime=rep(0,times=nInterim+1)
			isaDesignList[[i]]$dashboard$analysisTime[1]=time
			isaDesignList[[i]]$dashboard$currentSize=rep(0,times=nArms)
			isaDesignList[[i]]$dashboard$controlSize=0
		}
	}
	}
	return(isaDesignList)
}

# runSingleTrialSimulation
# Input: 
#	platformDesign: The building information for the platform design and ISAs
#	isaSimulationScenarios: The scenarios used to generate the response within ISAs
#	platformSimulationParameters: Computational parameters to guide the simulations
# Output:
#	isaDesignList: The final isaDesignList
#	fullData: the final data from the platform
#	allocationRates: optional allocation rates across ISAs for each time unit.

runSingleTrialSimulation=function(platformDesign,isaSimulationScenarios,platformSimulationParameters){

fullData=NULL
parkedData=NULL
subjData=NULL
isaSimulationParameterList=NULL
isaDesignList=getInitialDesignList(platformDesign)
dataSummary=initializeDataSummary(NULL,isaDesignList)
platformInactive=0
i=1
anyNew=TRUE
if(platformSimulationParameters$storeAllocRates){
	allocationRates=array(0,dim=c(length(platformDesign$entryTime),200))
}

k=0
kj=0

while(platformInactive<platformSimulationParameters$simulationSettings[1]){
	k=k+1
	if(anyNew){
		### Do all adjustments to isaDesignList only when something changes. Otherwise it is ok to work with isaDesignListTMP
		isaDesignList=prepareISADecisionResult(i,isaDesignList)
		if(platformSimulationParameters$useRandomizationList){
			isaDesignList=addISARandomizationList(isaDesignList,platformSimulationParameters,platformDesign)
		}
		
		nISA=length(isaDesignList)
		isaSimulationParameterList=updateSimulationParameters(isaDesignList,isaSimulationParameterList,isaSimulationScenarios,platformDesign)
		isaDesignListTMP=getFlatStructure(isaDesignList,isaSimulationParameterList)

		anyNew=FALSE
		updateAlloc=1
	}
	
	openEnrollmentPerISA=getVectorFromList(isaDesignListTMP,"activeISA",sumElements=TRUE)
	if(updateAlloc){
		platformDesign$allocationProbability=getAllocationProbability(isaDesignListTMP,platformDesign)
		if(platformSimulationParameters$useRandomizationList){
			nIsPerISA=sum(getVectorFromList(isaDesignListTMP,"currentSize",sumElements=TRUE))
			tmpList=getRandomizationListPlatform(platformDesign$allocationProbability,platformSimulationParameters)
			if(is.null(platformDesign$randomizationList)){
				platformDesign$randomizationList=tmpList
			} else {
				platformDesign$randomizationList=c(platformDesign$randomizationList[1:nIsPerISA],	
					getRandomizationListPlatform(platformDesign$allocationProbability,platformSimulationParameters))
			}
		}
		isaDesignListTMP=updateTargetSize(platformDesign$allocationProbability,isaDesignListTMP)
		updateAlloc=0
	}

	if(platformSimulationParameters$storeAllocRates){
		if((i%%201)==0){
			allocationRates=cbind(allocationRates,array(0,dim=c(length(platformDesign$entryTime),200)))
			allocationRates[(1:length(platformDesign$allocationProbability)),i]=platformDesign$allocationProbability
		} else {
			allocationRates[(1:length(platformDesign$allocationProbability)),i]=platformDesign$allocationProbability
		}
	}	else {
		allocationRates=NULL
	}
	
	if(any(openEnrollmentPerISA>0)){
				newPatientData=generateNewPatientData(i,isaDesignListTMP,platformDesign,platformSimulationParameters)
			if(is.null(parkedData)){
				parkedData=newPatientData
			} else {
				if(!is.null(newPatientData)){
					parkedData=rbind(parkedData,newPatientData)
				}
			}
			subjData=getSubjDataFile(newPatientData,subjData)
			subjDataTMP=getSubjDataFile(newPatientData,NULL)
		
			dataSummary=updateDataSummary(dataSummary,subjDataTMP)
			isaDesignListTMP=updateISAList(isaDesignListTMP,dataSummary)
			
			### Check if the number of open arms has changed. If so, the allocation ratio will need to be updated.
			openEnrollmentPerISA_New=getVectorFromList(isaDesignListTMP,"activeISA",sumElements=TRUE)
			if(sum(openEnrollmentPerISA)!=sum(openEnrollmentPerISA_New)){
				updateAlloc=1
			}
	} 

	# the full data set is not updated in every iteration to possibly save some time
	if(!is.null(parkedData)){
		if(dim(parkedData)[1]>platformSimulationParameters$simulationSettings[3]){
			if(is.null(fullData)){
				fullData=parkedData
			} else {
			fullData=rbindlist(list(fullData,parkedData))
			}
			parkedData=NULL
		}
	}
	if(k>=platformSimulationParameters$simulationSettings[2] & !is.null(subjData)){
		# only check all "k" days/weeks for new analyses
		isaDataTMP=getInterimDataPerISAFromSubjData(i,isaDesignListTMP,subjData)
		conductAnalysis=checkAnalysisMilestonesPlatform(isaDesignListTMP,i,isaDataTMP)
		if(any(conductAnalysis>0)){
			if(is.null(fullData)){
				fullData=parkedData 
			} else {
				fullData=rbindlist(list(fullData,parkedData))
			}
			parkedData=NULL
			isaData=getInterimDataPerISA(i,conductAnalysis,isaDesignListTMP,fullData)
			analysisResults=analysisWrapper(i,isaDesignListTMP,isaData,conductAnalysis)
			decisions=getDecisions(conductAnalysis,isaDesignListTMP,analysisResults)
			isaDesignListTMP=implementDecisions(i,conductAnalysis,isaDesignListTMP,analysisResults,decisions)

			isaDesignListTMP[[2]]=getControlAnalysis(isaDesignListTMP[[2]],fullData[which(fullData$isa==2 & fullData$visitTime<=i),])
			updateAlloc=1
		}
		k=0
	}
	addISA=addNewISA(i,isaDesignListTMP,isaDesignList,platformDesign)
	anyNew=addISA[[2]]>0
	if(anyNew){
		isaDesignList=addISA[[1]]
		dataSummary=initializeDataSummary(dataSummary,isaDesignList)
	}
	platformInactive=checkActivity(isaDesignListTMP,platformInactive)
	i=i+1
	if(i>platformSimulationParameters$simulationSettings[4]){
		print("NOT COMPLETE")
		break;
	}
}

return(list(isaDesignListTMP,fullData,allocationRates[,1:(i-1)]))
}

# checkActivity
# Input: 
#	isaDesignList: The isaDesignList
#	platformInactive: a counter how long the platform has been inactive
# Output:
#	integer counting how long no arms have been active
checkActivity=function(isaDesignList,platformInactive){
	nISA=length(isaDesignList)
	active=0
	for(i in seq(3,nISA)){
		if(isaDesignList[[i]]$activeISA){
			return(0)
		}
	}
	return((platformInactive+1))
}

