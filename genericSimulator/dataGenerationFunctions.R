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

# generateNewPatientData
# Input: 
#	time: Time in the platform
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# 	platformDesign: the overarching platform design elements.
# 	platformSimulationParameters:  list with simulation options.
# Output:
#	newSubjectData: data.table with complete data in long format for newly recruited subjects

# To-Do: Check format of newSubjectData (data.frame or data.table?)
generateNewPatientData=function(time,isaDesignList,platformDesign,platformSimulationParameters){
	nNewSubjects=recruitment(time,platformSimulationParameters$recruitmentPars,platformSimulationParameters$recruitmentType)
	## Generate a new data frame for the newly added subjects
	if(nNewSubjects>0){
		newSubjectData=getBaselineData(time,nNewSubjects,platformSimulationParameters$baselineSimulationPars)

	## Assign newly added subjects to ISA
	# ISA=0 => Longitudinal cohort
	# ISA=1 => Control cohort
	# ISA>1 => Investigative cohorts
	
	# New subject data is a data.table object
	
	### Get subjects to ISA
		newSubjectData=allocateToISA(newSubjectData,isaDesignList,platformDesign)

		#print(newSubjectData)
		newSubjectData=getAllSubjectsData(newSubjectData,isaDesignList,platformSimulationParameters,platformDesign)

		newSubjectData=newSubjectData[visitID<=dropOutTime]
		newSubjectData$entryTime=time
		newSubjectData$visitTime=time+newSubjectData$visitID

		if(dim(newSubjectData)[1]==0){
			return(NULL)
		} else {
			return(newSubjectData)
		}
	} else {
		return(NULL)
	}
}
	
# getBaselineData
# Input: 
#	time: Time in the platform
#	newSubjects: integer indicating the number of subjects newly recruited.
# 	parameters: parameters for defining covariates needs
#		nCovariates = the number of covariates
#		covariateType = the type of covariate
#		covariateNames = the name of the covariate (as used in the data)
#		parameter  = specific parameters defining the covariate distribution.
# Output:
#	outDat=array with one row per new patient and data on covariates in columns.
getBaselineData=function(t,newSubjects,parameters){
	outDat=cbind(as.numeric(paste(t,seq(1,newSubjects),sep=".")),t+runif(newSubjects))
	for(i in seq(1,parameters$nCovariates)){
		if(parameters$covariateType[i]=="Binary"){
			outDat=cbind(outDat,rbinom(newSubjects,1,parameters[[i+3]]))
		} else if(parameters$covariateType[i]=="Uniform"){
			outDat=cbind(outDat,runif(newSubjects,min=parameters[[i+3]][1],max=parameters[[i+3]][2]))
		} else if(parameters$covariateType[i]=="Normal"){
			outDat=cbind(outDat,rnorm(newSubjects,mean=parameters[[i+3]][1],sd=parameters[[i+3]][2]))
		} else if(parameters$covariateType[i]=="Site"){
			if(t<parameters[[i+3]][1]){
				outDat=cbind(outDat,rep(0,newSubjects))
			} else {
				outDat=cbind(outDat,rbinom(newSubjects,1,parameters[[i+3]][2]))
			}
		} else {
			outDat=cbind(outDat,getCustomBaselineCovariates(newSubjects,parameters[[i+3]]))
		}
	}
	colnames(outDat)=c("patID","patID_2",parameters$covariateNames)
	return(outDat)
}

# getAllSubjectsData
# Input: 
#	subjectData: array with data for subjects, including covariates.
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# 	platformSimulationParameters:  list with simulation options.
# 	platformDesign: the overarching platform design elements.
# Output:
#	outDat=data.table with the full newly generated subject data in long format.
getAllSubjectsData=function(subjectData,isaDesignList,platformSimulationParameters,platformDesign){
	nSubjects=dim(subjectData)[1]
	outMat=list()
	for(i in seq(1,nSubjects)){
		outMat=getSubjectData(subjectData[i,],isaDesignList,platformSimulationParameters,platformDesign)
		if(i==1){
			outDat=outMat
		} else {
			outDat=rbind(outDat,outMat)
		}
	}
	outDat=data.table(outDat)
	return(outDat)
}


# getSubjectArm
# Input: 
#	subjectData: a vector with baseline data for one specific subject.
#	isaDesign: a list containing a specific ISA designs to which the participant has been allocated.
# 	platformSimulationParameters:  list with simulation options.
# 	platformDesign: the overarching platform design elements.
# Output:
#	arm: the arm within the ISA to which the participant has been randomized.
getSubjectArm=function(subjectData,isaDesign,platformSimulationParameters,platformDesign){
	subjectDataIn=subjectData

	### ISA assignment is stored in last element of subjectDataIn.
	currentISA=subjectDataIn[length(subjectDataIn)-1]
	if(currentISA==2){
		arm=1
		return(arm)
	}
	if(platformDesign$randomizer=="List"){
		### Allocate according to randomization list, if available.
		arm=isaDesign$randomizationList[subjecData$pat_ISA_ID]
	} else if(platformDesign$randomizer=="RAR"){
		### Allocate proportional to the targeted allocation Ratio.
		### Tuning parameter: platformDesign$randomizationParameter
		### Correct in the end for too low randomization probability
		### TBC: Should this be done earlier when calculating the allocation ratios?
		if(isaDesign$allocationRule=="constant"){
			if(length(isaDesign$targetSize)==1){
				return(1)
			}
			targetedRandomizationRate=rep(1/length(isaDesign$allocProb),length(isaDesign$allocProb))		
			openArms=which(isaDesign$currentSize<isaDesign$targetSize)		
			if(any(isaDesign$currentSize>=isaDesign$targetSize)){
				closedArms=which(isaDesign$currentSize>=isaDesign$targetSize)
				targetedRandomizationRate[closedArms]=0
			}
			
			if(any(targetedRandomizationRate[openArms]<platformDesign$minProb)){
					indexes=which(targetedRandomizationRate[openArms]<platformDesign$minProb)
					cut=sum(targetedRandomizationRate[openArms]<platformDesign$minProb)*platformDesign$minProb
					targetedRandomizationRate[indexes]=platformDesign$minProb
					targetedRandomizationRate[-indexes]=(1-cut)*targetedRandomizationRate[-indexes]
			} 

			if(sum(isaDesign$currentSize)==0){
				currentRandomizationRate=rep(1/length(isaDesign$currentSize),times=length(isaDesign$currentSize))
			} else {
				currentRandomizationRate=isaDesign$currentSize/sum(isaDesign$currentSize)
			}
			
			arm=randomizeSubject(currentRandomizationRate,
							targetedRandomizationRate,platformSimulationParameters$randomizationParameter)
			
		} else {			
			openArms=which(isaDesign$currentSize<isaDesign$targetSize)
			allocProb=isaDesign$allocProb
			
			if(any(isaDesign$currentSize>=isaDesign$targetSize)){
				closedArms=which(isaDesign$currentSize>=isaDesign$targetSize)
				allocProb[closedArms]=0
			}
			
			if(sum(isaDesign$allocProb[openArms])==0){
				isaDesign$allocProb[openArms]=1
			}
			
			targetedRandomizationRate=sqrt(allocProb)/sum(sqrt(allocProb))
			randProb=cumsum(targetedRandomizationRate/sum(targetedRandomizationRate))
			tmp=runif(1)
			arm=min(which(randProb>tmp))
		}		
	} else {
		### Allocate within ISA according to to targeted final randomization probability
		if(sum(isaDesign$currentSize==0)){
			currentRandomizationRate=rep(1/length(isaDesign$currentSize),times=length(isaDesign$currentSize))
		} else {
			currentRandomizationRate=isaDesign$currentSize/sum(isaDesign$currentSize)
		}
		targetedRandomizationRate=isaDesign$targetSize/sum(isaDesign$targetSize)
		targetedRandomizationRate[which(isaDesign$targetSize<=isaDesign$currentSize)]=0
		targetedRandomizationRate=targetedRandomizationRate/sum(targetedRandomizationRate)
		
		arm=randomizeSubject(currentRandomizationRate,
						targetedRandomizationRate,platformSimulationParameters$randomizationParameter)	
	}
	return(arm)
}

# getSubjectData
# Input: 
#	subjectData: a vector with baseline data for one specific subject.
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# 	platformSimulationParameters:  list with simulation options.
# 	platformDesign: the overarching platform design elements.
# Output:
#	outMat: an array containing the complete simulated data for a subject in long format.
getSubjectData=function(subjectData,isaDesignList,platformSimulationParameters,platformDesign){
	subjectDataIn=subjectData
	currentISA=subjectDataIn[length(subjectDataIn)-1]
	isaDesign=isaDesignList[[currentISA]]
	
	arm=getSubjectArm(subjectData,isaDesign,platformSimulationParameters,platformDesign)
	
	visitVec=isaDesign$visitVec
	
	if(any(platformSimulationParameters$endpoint=="customized")){
		tmpDat=customizedData(arm,visitVec,isaDesign$endpointParameters,isaDesign)
		dropOutTime=getDropOutTime(isaDesign$dropOutEventRate[arm],tmpDat)
		# First column hard coded: Must be visitID
		censoring=tmpDat[,1]>=dropOutTime	
		nRows=dim(tmpDat)[1]
		endpointNames=colnames(tmpDat)
		
		tmpDat=cbind(tmpDat,rep(dropOutTime,times=nRows),censoring,rep(arm,times=nRows),rep(isaDesign$trtSpecifics[arm],times=nRows))
		
		tmp=as.numeric(subjectDataIn)
		outMat=repRows(tmp,nRows)
		outMat=cbind(outMat,tmpDat)
	
		colnames(outMat)=c(names(subjectDataIn),endpointNames,"dropOutTime","censored","arm","trtSpecifics")
	} else {
	
		for(j in seq(1,platformSimulationParameters$nEndpoints)){
			if(j==1){
				tmpDat=getEndpointData(j,arm,visitVec,platformSimulationParameters$endpoint[j],isaDesign$endpointParameters[[j]])
			} else {
				newDat=getEndpointData(j,arm,visitVec,platformSimulationParameters$endpoint[j],isaDesign$endpointParameters[[j]])
				tmpDat=cbind(tmpDat,newDat)
			}
		}
	
	endpointNames=colnames(tmpDat)
	dropOutTime=getDropOutTime(isaDesign$dropOutEventRate[arm])
	censoring=visitVec>=dropOutTime

	nRows=dim(tmpDat)[1]
	tmpDat=cbind(tmpDat,rep(dropOutTime,times=nRows),censoring,visitVec,rep(arm,times=nRows),rep(isaDesign$trtSpecifics[arm],times=nRows))
	
	tmp=as.numeric(subjectDataIn)
	outMat=repRows(tmp,nRows)
	outMat=cbind(outMat,tmpDat)
	colnames(outMat)=c(names(subjectDataIn),endpointNames,"dropOutTime","censored","visitID","arm","trtSpecifics")
	}
	return(outMat)
}

# getDropOutTime
# Input: 
#	parameters = yearly drop-out rate (assuming data simulated weekly!)
# Output:
#	out = simulated drop-out time.

# To-Do: Ensure that drop-out time can be customized to time unit used for simulations.
getDropOutTime=function(parameters,tmpDat){
	rate=-log(1-parameters)/52
	out=rexp(1,rate=rate)
}

# getEndpointData
# Input: 
#	index: indicating the "number of the endpoint"
#	arm:   arm within ISA 
#	visitVec: assessment schedule
#	endpoint: type of endpoint
#	parameters: parameters used to define distribution of endpoint
#	covariates: covariates from patient to be taken into account in simulation
# Output:
#	out = array with one row per visit.

# To-Do: Test all endpoint types
getEndpointData=function(index,arm,visitVec,endpoint,parameters,covariates){
	if(endpoint=="cont"){
		out=array(0,dim=c(length(visitVec),1))
		indParVec=rmvnorm(1,mean=parameters$meanPars,sigma=parameters$D)
		out=t(t(indParVec[1]+visitVec*indParVec[2]+rnorm(length(visitVec))*parameters$sd))
		colnames(out)="contResp"
	} else if(endpoint=="contBinary"){
		out=array(0,dim=c(length(visitVec),3))
		indParVec=rmvnorm(1,mean=parameters$meanPars,sigma=parameters$D)
		out[,1]=indParVec[1]+visitVec*indParVec[2]+rnorm(length(visitVec))*parameters$sd
		out[,2]=out[,1]<parameters$cutoff
		if(any(out[,2]==1)){
			out[min(which(out[,2]==1)):length(visitVec),3]=1
		}
		colnames(out)=c("contResp","binResp","binResp2")
	} else if(endpoint=="timeToEvent"){
		out=array(0,dim=c(length(visitVec),2))
		lambdaPar=exp(log(-log(1-parameters$eventRate[arm]))/parameters$shapeParameter)/parameters$followUp
		eventTime=rweibull(1,shape=parameters$shapeParameter,scale=1/lambdaPar)
		out[,1]=rep(eventTime,length(visitVec))
		out[,2]=0
		out[which(visitVec>=eventTime),2]=1
		colnames(out)=c("eventTime","binResp")
	} else if(endpoint=="binary"){
		
		out=array(0,dim=c(length(visitVec),1))
		if(length(visitVec)>2){
			logits=rmvnorm(1,mean=parameters$meanPars,sigma=parameters$D)
		} else {
			logits=c(-9999,parameters$meanPars+parameters$D*rnorm(1))
		}
		out[,1]=rbinom(length(visitVec),1,prob=exp(logits)/(1+exp(logits)))
		colnames(out)="binResp"
	} else if(endpoint=="weibull"){
		out=array(0,dim=c(length(visitVec),2))
		lambdaPar1=exp(log(-log(1-parameters$eventRate1[arm]))/parameters$shapeParameter1[arm])/parameters$followUp1[arm]
		eventTime=rweibull(1,shape=parameters$shapeParameter1[arm],scale=1/lambdaPar1)+parameters$shift[arm]
			
		out[,1]=rep(eventTime,length(visitVec))
		out[,2]=0
		out[which(visitVec>=eventTime),2]=1
		colnames(out)=paste(c("eventTime","binResp"),index,sep=".")
	} else if(endpoint=="piecewiseExponential"){
		out=array(0,dim=c(length(visitVec),2))
		lambdaPar1=exp(log(-log(1-parameters$eventRate1[arm]))/parameters$shapeParameter1)/parameters$followUp1
		lambdaPar2=exp(log(-log(1-parameters$eventRate2[arm]))/parameters$shapeParameter2)/parameters$followUp2
		eventTime=rweibull(1,shape=parameters$shapeParameter1,scale=1/lambdaPar1)
		
		if(eventTime>parameters$duration[arm]){
			eventTime=parameters$duration[arm]+rweibull(1,shape=parameters$shapeParameter2,scale=1/lambdaPar2)
		}
		out[,1]=rep(eventTime,length(visitVec))
		out[,2]=0
		out[which(visitVec>=eventTime),2]=1
		colnames(out)=c("eventTime","binResp")
	} else if(endpoint=="piecewiseExponential2"){
		out=array(0,dim=c(length(visitVec),2))
		t=parameters$followUp
		lambdaPar1=-log(1-parameters$eventRate[1])/t[1]
		lambdaPar2=-(log(1-parameters$eventRate[2])+lambdaPar1*t[1])/(t[2]-t[1])
		lambdaPar3=-(log(1-parameters$eventRate[3])+lambdaPar1*t[1]+lambdaPar2*(t[2]-t[1]))/(t[3]-t[2])
		
		eventTime=rexp(1,rate=lambdaPar1)
		if(eventTime>t[1]){
			eventTime=t[1]+rexp(1,rate=lambdaPar2)
		} 
		if(eventTime>t[2]){
			eventTime=t[2]+rexp(1,rate=lambdaPar3)
		}
		
		out[,1]=rep(eventTime,length(visitVec))
		out[,2]=0
		out[which(visitVec>=eventTime),2]=1
		colnames(out)=c("eventTime","binResp")
	} else if(endpoint=="binaryTimeToEvent"){
		out=array(0,dim=c(length(visitVec),3))
		responder=rbinom(1,size=1,prob=parameters$responseRate[arm])
		if(responder==1){
			lambdaPar1=exp(log(-log(1-parameters$eventRate[arm]))/parameters$shapeParameter[arm])/parameters$followUp[arm]
			eventTime=rweibull(1,shape=parameters$shapeParameter[arm],scale=1/lambdaPar1)+parameters$shift[arm]
		} else {
			eventTime=Inf
		}
		out[,1]=rep(eventTime,length(visitVec))
		out[,2]=0
		out[which(visitVec>=eventTime),2]=1
		if(any(out[,2]==1)){
			out[,3]=visitVec[min(which(visitVec>=eventTime))]
		} else {
			out[,3]=Inf
		}
		colnames(out)=paste(c("eventTime","binResp","censorTime"),index,sep=".")
	} else {
		return(rep(NA,length(visitVec)))
	}
	return(out)
}

