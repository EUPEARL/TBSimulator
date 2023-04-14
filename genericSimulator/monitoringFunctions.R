##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions which check whether an interim analysis should be conducted.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# tiggerByEvents
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on entryTime, event time and drop out time
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered, if the required number of total events has been reached (not control group only)
tiggerByEvents=function(isaDesign,t,data){
	conductAnalysis=0
	eventTime=data$eventTime
	tmp1=(eventTime<data$dropOutTime)
	tmp2=(eventTime+data$entryTime<t)
	nEvents=sum(tmp1*tmp2)
	analysisCondition=nEvents>=isaDesign$interimCutoff
	if(any(analysisCondition & isaDesign$interimPassed==0)){
		conductAnalysis=max(which(analysisCondition & isaDesign$interimPassed==0))
	}
	return(conductAnalysis)
}

# triggerByInformationCompleted
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered, if the required relative information (completers) has been reached. Information defined as % of active subjects.
triggerByInformationCompleted=function(isaDesign,t,data){
	conductAnalysis=0
	subData=data[which(data$followUp>=max(isaDesign$visitVec)),]
	subData=subData[which(subData$arm!=0),]
	relativeInformation=dim(subData)[1]/sum(isaDesign$targetSize)
	analysisCondition=relativeInformation>=isaDesign$interimCutoff
	if(any(analysisCondition & isaDesign$interimPassed==0)){
		conductAnalysis=max(which(analysisCondition & isaDesign$interimPassed==0))
	} 
	return(conductAnalysis)
}

# triggerByTime
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered, if the ISA has been open for long enough. Interim analysis are conducted by "cutoff"
triggerByTime=function(isaDesign,t,data){
	conductAnalysis=0
	analysisCondition=(t-isaDesign$startTime)>=isaDesign$interimCutoff
	if(any(analysisCondition & isaDesign$interimPassed==0)){
		conductAnalysis=max(which(analysisCondition & isaDesign$interimPassed==0))
	} 
	return(conductAnalysis)
}

# triggerByInterval
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered repeatedly every "interimCutoffs"-weeks
triggerByInterval=function(isaDesign,t,data){
	conductAnalysis=0
	startTime=isaDesign$analysisTime[1]
	burnIn=(t-startTime)
	analysisCondition=(t-startTime)%%isaDesign$interimCutoff==0
	openEnrollment=isaDesign$activeISA

	if(analysisCondition & openEnrollment){
		conductAnalysis=isaDesign$stage
	} else if(!openEnrollment){
		lpi=max(data[which(data$arm!=0),]$entryTime)
		if(lpi-isaDesign$followUp>0){
			conductAnalysis=isaDesign$stage
		}
	}
	return(conductAnalysis)
}

# triggerByInformationEnrolled
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered based on relative number of targeted subjects enrolled on active arms.
triggerByInformationEnrolled=function(isaDesign,t,data){
	conductAnalysis=0
	subData=data[which(data$arm!=0),]
	relativeInformation=dim(subData)[1]/sum(isaDesign$targetSize)
	analysisCondition=relativeInformation>=isaDesign$interimCutoff
	if(any(analysisCondition & isaDesign$interimPassed==0)){
		conductAnalysis=max(which(analysisCondition & isaDesign$interimPassed==0))
	} 
	return(conductAnalysis)
}


# triggerByFollowUp
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered based on number of active participants who reached a certain targeted followUp.
triggerByFollowUp=function(isaDesign,t,data){
	conductAnalysis=0
	subData=data[which(data$arm!=0),]$followUp
	
	nStages=length(isaDesign$interimCutoff)
	analysisCondition=rep(0,nStages)
	
	for(i in seq(1,nStages)){
		followUp=isaDesign$interimFollowUp[i]
		analysisCondition[i]=sum(subData>=followUp)>=isaDesign$interimCutoff[i]
	}

	if(any(analysisCondition & isaDesign$interimPassed==0)){
		conductAnalysis=max(which(analysisCondition & isaDesign$interimPassed==0))
	} 
	return(conductAnalysis)
}

# checkAnalysisMilestones
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	interim analysis, which is to be conducted or "0", if no interim needed.
# 	interim triggered based on one of multiple triggers. If the ISA is already complete, the value "0" will be provided.
# 	If ISA is inactive and not complete, interim analysis may still be triggered by events, followUp or custom trigger.
#	  For other approaches, wait for all remaining data to complete and conduct the final analysis.
#	Note: If customTrigger used, will need to specify this customized function in "customFunctions.R".
checkAnalysisMilestones=function(isaDesign,t,data){
	if(length(isaDesign$completed)>1){
	break;
	}
	
	if(isaDesign$completed){
		return(0)
	}
	
	### Case: not enrolling, but final analysis not reached yet
	if(!isaDesign$activeISA & !isaDesign$completed){
		if(isaDesign$interimTrigger=="events"){
			return(tiggerByEvents(isaDesign,t,data))
		}  else if(isaDesign$interimTrigger=="followUp"){
			return(triggerByFollowUp(isaDesign,t,data))	
		} else if(isaDesign$interimTrigger=="customTrigger"){
			return(customTrigger(isaDesign,t,data))
		} else {
			conductAnalysis=0
			subData=data[which(data$arm!=0),]
			if(t>max(subData$entryTime)+max(isaDesign$visitVec)){
				conductAnalysis=1
			}
			return(conductAnalysis)
		}
	} else {
			if(isaDesign$interimTrigger=="informationCompleted"){
				return(triggerByInformationCompleted(isaDesign,t,data))
			} else if(isaDesign$interimTrigger=="informationEnrolled"){
				return(triggerByInformationEnrolled(isaDesign,t,data))
			} else if(isaDesign$interimTrigger=="time"){
				return(triggerByTime(isaDesign,t,data))	
			} else if(isaDesign$interimTrigger=="perpetual"){
				return(triggerByInterval(isaDesign,t,data))	
			} else if(isaDesign$interimTrigger=="followUp"){
				return(triggerByFollowUp(isaDesign,t,data))	
			} else if(isaDesign$interimTrigger=="events"){
				return(tiggerByEvents(isaDesign,t,data))
			} else if(isaDesign$interimTrigger=="customTrigger"){
				return(customTrigger(isaDesign,t,data))
			}
	}
}

# checkAnalysisMilestonesPlatform
# Input: 
#	isaDesign: isaDesignList, containing informations on cutoffs and stage
#	t: Time in platform
#	data: containing information on followUpTime
# Output:
#	Wrapper going through ISAs and checking for need of doing interim analysis. Returns vector describing need for IA.
checkAnalysisMilestonesPlatform=function(isaDesignList,t,isaData){
		nISA=length(isaDesignList)
		conductAnalysis=c()
		k=1
		for(i in seq(3,nISA)){
			isaDataTMP=isaData[[k]]
			if(dim(isaDataTMP)[1]>0){
				conductAnalysis[k]=checkAnalysisMilestones(isaDesignList[[i]],t,isaData[[k]])
			} else {
				conductAnalysis[k]=0
			}
			k=k+1
		}
	return(conductAnalysis)
}

# initializeDataSummary
# Input: 
#	dataSummary: The dataSummary-Object is a table, which summarizes data on the sample size enrolled. It has 5 rows: ISA,ARM within ISA,minimum entry time, maximum entry time, sample size.
#	isaDesignList: Contains information on the number of ISAs, arms and targeted size
# Output:
#	dataSummary
initializeDataSummary=function(dataSummary,isaDesignList){
	nISA=length(isaDesignList)
	nArms=0
	for(i in seq(1,nISA)){
		nArms=nArms+length(isaDesignList[[i]]$designElements$targetSize)
	}
	if(is.null(dataSummary)){
		dataSummary=array(0,dim=c(5,nArms))
	} else {
		oldArms=dim(dataSummary)[2]
		addedColumns=nArms-oldArms
		dataSummary=cbind(dataSummary,array(0,dim=c(5,addedColumns)))
	}
	
	j=1
	for(i in seq(1,nISA)){
		newArms=length(isaDesignList[[i]]$designElements$targetSize)
		if(dataSummary[1,j]==0){
			dataSummary[1,j:(j+newArms-1)]=i
			if(i<3){
				dataSummary[2,j]=0
			} else {
				dataSummary[2,j:(j+newArms-1)]=seq(1,newArms)
			}
		}
		j=j+newArms
	}
	return(dataSummary)
}

# updateDataSummary
# Input: 
#	dataSummary: The dataSummary-Object is a table, which summarizes data on the sample size enrolled. It has 5 rows: ISA,ARM within ISA,minimum entry time, maximum entry time, sample size.
#	newData: if new data has been generated, contains new patient level data with information entryTime, arm and isa.
# Output:
#	updated dataSummary (5 rows: ISA,ARM within ISA,minimum entry time, maximum entry time, sample size, one column for each ISA and arm)
updateDataSummary=function(dataSummary,newData){
	if(is.null(newData)){
		return(dataSummary)
	}
	
	uniqueISA=as.numeric(unique(newData$isa))

	for(i in seq(1,length(uniqueISA))){
		columns=which(dataSummary[1,]==uniqueISA[i])
		subDat=newData[which(newData$isa==uniqueISA[i]),]
		entry=min(subDat$entryTime)
		arms=unique(subDat$arm)
		for(j in seq(1,length(arms))){
			if(min(dataSummary[3,columns[arms[j]]])!=0){
				dataSummary[3,columns[arms[j]]]=min(entry,min(dataSummary[3,columns[arms[j]]]))
			} else {
				dataSummary[3,columns[arms[j]]]=min(entry)
			}
			dataSummary[4,columns[arms[j]]]=max(entry,max(dataSummary[4,columns[arms[j]]]))
			dataSummary[5,columns[arms[j]]]=dataSummary[5,columns[arms[j]]]+length(which(subDat$arm==arms[j]))
		}
	}
	return(dataSummary)
}


