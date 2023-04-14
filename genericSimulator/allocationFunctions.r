##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions to support allocation of patients across ISAs.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# checkAdmissibleISAs
# Input: 
#	patientDat: vector with data for the patient. In particular a population identifier. If this population is part of the ISA, the patient can be included.
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# Output:
#	out: vector indicating whether assignment of patient to isa is appropriate.
checkAdmissibleISAs=function(patientDat,isaDesignList){
	### To be filled with rules for different inclusion/exclusion criteria
	#return(rep(1,length(isaDesignList)))
	
	out=rep(0,length(isaDesignList))
	for(i in seq(1,length(isaDesignList))){
		if(i==1){
			out[i]=1e-5
		} else {
			if(is.null(isaDesignList[[i]]$subPopulation)){
				out[i]=1
			} else {
				column=which(names(patientDat)==getEnrichmentVariable())
				out[i]=any(isaDesignList[[i]]$subPopulation==patientDat[column])
			} 
		} 
	}
	return(out)
}

# getISAPatientIndex
# Input: 
#	isaVec: new assignments to the ISAs
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# Output:
#	outVec: counter indicating for each added patient the sequence in the assigned isa
getISAPatientIndex=function(isaVec,isaDesignList){
	alloc=unique(isaVec)
	outVec=rep(0,length(isaVec))
	#print(isaVec)
	for(i in alloc){
		nToISA=length(which(isaVec==i))
		currentSize=sum(isaDesignList[[i]]$currentSize)
		#print(c("currentSize",currentSize))
		outVec[which(isaVec==i)]=seq((currentSize+1),(currentSize+nToISA))
	}
	return(outVec)
}

# allocateToISA
# Input: 
#	newPatientData: An array with one row of data for each newly added patient
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# 	platformDesign: the overarching platform design elements.
# Output:
#	newPatientData: The input array newPatientData with additional information on the ISA assignment and the patientID within the ISA.

# ToDo: Need to check and optimize function for use of "randomizationList" and two-stage randomization.
allocateToISA=function(newPatientData,isaDesignList,platformDesign){
	nSubjects=dim(newPatientData)[1]
	nTargetPerISA=getVectorFromList(isaDesignList,"targetSize",sumElements=TRUE)
	nIsPerISA=getVectorFromList(isaDesignList,"currentSize",sumElements=TRUE)
	openEnrollmentPerISA=getVectorFromList(isaDesignList,"activeISA",sumElements=TRUE)
	priorAllocationProbability=platformDesign$allocationProbability/sum(platformDesign$allocationProbability)

	isaVec=c()
	
	if(any(names(platformDesign)=="randomizationList")){
		for(i in seq(1,nSubjects)){
			isaVec[i]=platformDesign$randomizationList[sum(nIsPerISA)+i]
		}
	} else {
		for(i in seq(1,nSubjects)){			
			indexVec=checkAdmissibleISAs(newPatientData[i,],isaDesignList)
			allocationProbability=cumsum((openEnrollmentPerISA*priorAllocationProbability*indexVec)/sum((openEnrollmentPerISA*priorAllocationProbability*indexVec)))
			
			tmp=runif(1)

			if(!any(openEnrollmentPerISA*priorAllocationProbability*indexVec>0)){
				isaVec[i]=0
			} else {
				isaVec[i]=min(which(allocationProbability>tmp))
			}
		}
	}
	
	colnamesTMP=colnames(newPatientData)
	isaPatIndex=getISAPatientIndex(isaVec,isaDesignList)
	newPatientData=cbind(newPatientData,isaVec,isaPatIndex)	
	colnames(newPatientData)=c(colnamesTMP,"isa","pat_ISA_ID")
	return(newPatientData)
}

# getAllocationProbability
# Input: 
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
# 	platformDesign: the overarching platform design elements.
# Output:
#	allocProb: Vector with randomization probability to each ISA.
getAllocationProbability=function(isaDesignList,platformDesign){
	if(platformDesign$randomizer=="RAR"){
		return(getCustomizedAllocationProbability(isaDesignList,platformDesign))
	} else {
		nISA=length(isaDesignList)
		for(i in seq(1,nISA)){
			openEnrollmentPerISA=getVectorFromList(isaDesignList,"activeISA",sumElements=TRUE)
		}
		nActiveArms=c(0,0)
		nArms=c(0,0)
		nActiveISA=c(0,0)
		for(i in seq(3,nISA)){
			if(openEnrollmentPerISA[i]){
				activeArms=isaDesignList[[i]]$currentSize<isaDesignList[[i]]$targetSize
				nActiveArms[i]=sum(activeArms)
				nArms[i]=length(activeArms)
				nActiveISA[i]=1
			} else {
				nActiveArms[i]=0
				nArms[i]=0
				nActiveISA[i]=0
			}
		}
		if(platformDesign$allocationRule[1]=="perISA"){
			if(platformDesign$allocationRule[2]=="SR"){
				rate=c(0,sqrt(sum(nActiveISA)),nActiveISA[-c(1:2)])
			} else {
				rate=c(0,as.numeric(platformDesign$allocationRule[4]),nActiveISA[-c(1:2)])
			}
		} else if(platformDesign$allocationRule[1]=="perArm"){
			if(platformDesign$allocationRule[2]=="SR"){
				rate=c(0,sqrt(sum(nArms)),nArms[-c(1:2)])
			} else {
				rate=c(0,as.numeric(platformDesign$allocationRule[4]),nArms[-c(1:2)])
			}
		} else if(platformDesign$allocationRule[1]=="perActiveArm"){
			if(platformDesign$allocationRule[2]=="SR"){
				rate=c(0,sqrt(sum(nActiveArms)),nActiveArms[-c(1:2)])
			} else {
				rate=c(0,as.numeric(platformDesign$allocationRule[4]),nActiveArms[-c(1:2)])
			}
		}
		if(sum(rate)==0){
			#print("Warning: No test intervention, assigning to control")
			return(c(0,1,rate[-c(1:2)]))
		}
		rate=rate
		ctrlAlloc=rate[2]/sum(rate)
		ctrlBoundary=as.numeric(platformDesign$allocationRule[3])
		if(ctrlAlloc<ctrlBoundary){
			ctrlAlloc=ctrlBoundary
		}
		
		if(sum(rate[-c(1:2)])!=0){	
			allocProb=rate[-c(1:2)]/sum(rate[-c(1:2)])
			allocProb=c(0,ctrlAlloc,(1-ctrlAlloc)*allocProb)
		} else {
			allocProb=c(0,1,rate[-c(1:2)])
		}
		return(allocProb)
	}
}