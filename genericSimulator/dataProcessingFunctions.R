##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions which process data for monitoring and use in interim analyses
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# getInterimDataPerISA
# Input: 
#	time: Time in the platform
#	conductAnalysis: a vector indicating which analysis is to be conducted per ISA. If "0", no analysis will be conducted.
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
#	fullData:	the full data frame in long format (one row per visit per patient)
# Output:
#	isaData: a list with compiled data for each ISA with interim analysis
getInterimDataPerISA=function(time,conductAnalysis,isaDesignList,fullData){
	nISA=length(isaDesignList)
	controlData=fullData[fullData$isa==2&fullData$visitTime<=time,]
	if(dim(controlData)[1]>0){
		controlData$arm=0
	}
	controlDataEntryTime=controlData$entryTime
	isaData=list()
	k=1
	for(i in seq(3,nISA)){
		if(conductAnalysis[[k]]){
			indexes=which(fullData$isa==i&fullData$visitTime<=time)
			isaData[[k]]=fullData[indexes,]
			if(dim(isaData[[k]])[1]>0){
				startDate=min(isaData[[k]]$entryTime)
				endDate=max(isaData[[k]]$entryTime)+3
				controlDataTMP=getControlData(controlData,startDate,endDate,isaDesignList[[i]],controlDataEntryTime)
				isaData[[k]]=rbind(isaData[[k]],controlDataTMP)
			}
		} else {
			isaData[[k]]=NULL
		}
			k=k+1
	}
	return(isaData)
}

# getInterimDataPerISAFromSubjData
# Input: 
#	time: Time in the platform
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
#	subjData:	the full data across ISAs, but in a reduced format with just one row per patient
# Output:
#	isaData: a list with reduced data file for each ISA with one row per patient 
getInterimDataPerISAFromSubjData=function(time,isaDesignList,subjData){

	nISA=length(isaDesignList)
	subjData$followUp=time-subjData$entryTime
	
	controlData=subjData[subjData$isa==2,]
	if(dim(controlData)[1]>0){
		controlData$arm=0
	}
	controlDataEntryTime=controlData$entryTime
	isaData=list()
	k=1
	for(i in seq(3,nISA)){
				indexes=which(subjData$isa==i)
				isaData[[k]]=subjData[indexes,]

				if(dim(isaData[[k]])[1]>0){
					startDate=min(isaData[[k]]$entryTime)
					endDate=max(isaData[[k]]$entryTime)+3
					controlDataTMP=getControlData(controlData,startDate,endDate,isaDesignList[[i]],controlDataEntryTime)
					isaData[[k]]=rbind(isaData[[k]],controlDataTMP)
				}
			k=k+1
	}
	return(isaData)
}

# getControlData
# Input: 
#	controlData: the full control data
#	startDate: 	FPI for an ISA
#	endData:	LPI for an ISA + 3 weeks
#	isaDesign:  one specific ISA design
#	entryTime:	vector of entry times of the control group.
# Output:
#	controlData: Selected control data for the ISA of interest (concurrent vs. non-concurrent, vs. subpopulations of interest)

# To-Do: 
##	Need to store the test than pool decision properly in the isaDesign
##	Need to adjust for proper implementation of subgroups, likely via API.
getControlData=function(controlData,startDate,endDate,isaDesign,entryTime){
	if(isaDesign$controlType=="concurrentControl"){
		indexes=(entryTime>=startDate & entryTime<endDate)
		controlData=(controlData[indexes,])
		if(dim(controlData)[1]==0){
			return(controlData)
		}
	} else if(isaDesign$controlType=="allControl"){
		controlData=controlData
	} else if(isaDesign$controlType=="testThanPool"){
		if(isaDesign$decisions[1,6]==1){
			controlData=controlData	
		} else {
			indexes=(entryTime>=startDate & entryTime<endDate)
			controlData=(controlData[indexes,])
			if(dim(controlData)[1]==0){
				return(controlData)
			}
		}
	}
	if(!is.null(isaDesign$subPopulation)){
		subpopulation=controlData$Site
		out=c()
		for(i in seq(1,length(subpopulation))){
			if(any(isaDesign$subPopulation==subpopulation[i])){
				out[i]=1
			} else {
				out[i]=0
			}
		}
		controlData=controlData[which(out==1),]
	}
	return(controlData)
}


# getSubjDataFile
# Input: 
#	patientData: new patient data, which is to be added to the subject data file. The last visit per patient is stored in the subjData file.
#		at minimum stored: patID_2, isa, arm, visitTime, entryTime
#		additional: eventTimes, responder status, dropOutTime, event censor enrichment variables
#	subjData: the subject data file with one row per patient.
# Output:
#	subjData: the updated subject data file with one row per patient.

# To-Do: 
##	Make the columns to be added to subjData file generic via API.
getSubjDataFile=function(patientData,subjData){
	if(is.null(patientData)){
		return(subjData)
	}
	columns=c("patID_2","isa","arm","visitTime","entryTime")
	if(any(names(patientData)=="eventTime")){
		columns=c(columns,"eventTime","binResp","dropOutTime")
	}	
	if(any(names(patientData)=="eventTime.1")){
		columns=c(columns,"eventTime.1","binResp.1","dropOutTime")
	}				
	if(any(names(patientData)=="eventTime.2")){
		columns=c(columns,"eventTime.2","binResp.2","dropOutTime")
	}				
	if(any(names(patientData)=="TTCC")){
		columns=c(columns,"TTCC","TTCC.1","TTCC.Censor","event","dropOutTime")
	}				
	if(any(names(patientData)==getEnrichmentVariable())){
		columns=c(columns,getEnrichmentVariable())
	}				
	patID=patientData$patID_2
	indexVec=getLastVisitIndexes(patID)
	dat=patientData[indexVec,columns,with=FALSE]
	
	if(is.null(subjData)){
		subjData=dat
	} else {
		subjData=rbind(subjData,dat)
	}
	return(subjData)
}