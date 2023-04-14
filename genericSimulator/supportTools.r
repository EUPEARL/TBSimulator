##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions supporting various small operations within the platform simulator
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# allocationFunction
# Input: 
#	pIs: Current randomization ratio
#	pTarget: Targeted randomization ratio
# 	gamma: Parameter tuning randomization (Small => Fixed sequence. Large => Random
# Output:
#	gOut= a weight for the respective arm
allocationFunction=function(pIs,pTarget,gamma){
	if(pIs==0){
		#return(1-1e-5)
		return(pTarget)
	} else if(pIs==1){
		return(1e-5)
	} else {
		tmp=(pTarget/pIs)^gamma
		gOut=(pTarget*tmp)/(pTarget*tmp+(1-pTarget)*((1-pTarget)/(1-pIs))^gamma)
		return(gOut)
	}
}

# randomizeMultiArm
# Input: 
#	pIs: Current randomization ratio
#	pTarget: Targeted randomization ratio
# 	gamma: Parameter tuning randomization (Small => Fixed sequence. Large => Random)
# Output:
#	arm= The final arm, which is randomized to
randomizeMultiArm=function(pIs,pTarget,gamma){
	nArms=length(pIs)
	if(nArms==1){
		return(1)
	} else {
		allocFunc=c()
		for(i in seq(1,nArms)){
			allocFunc[i]=allocationFunction(pIs[i],pTarget[i],gamma)
		} 
		randProb=cumsum(allocFunc/sum(allocFunc))
		tmp=runif(1)
		armOut=min(which(randProb>tmp))
		return(armOut)
	}
}
	
# randomizeSubject
# Input: 
#	currentRate: Current randomization ratio
#	targetRate: Targeted randomization ratio
# 	randomizationParameter: Parameter tuning randomization (Small => Fixed sequence. Large => Random)
# Output:
#	armOut= The final arm, which is randomized to
randomizeSubject=function(currentRate,targetRate,randomizationParameter){
	armOut=randomizeMultiArm(currentRate,targetRate,randomizationParameter)
}


# getRandomizationListPlatform
# Input: 
#	allocationProbability: vector with allocation probability to each arm
#	simulationParameters: list containing parameters used to define randomization list
#		blockSize
# Output:
# 	Vector of assignments
getRandomizationListPlatform=function(allocationProbability,simulationParameters){
  ratio=ceiling(allocationProbability*simulationParameters$blockSize)
  adjBlockSize=sum(ratio)
  groupNames=as.character(seq(1,length(ratio)))
  if(any(ratio==0)){
    groupNames=groupNames[which(ratio!=0)]
    ratio=ratio[which(ratio!=0)]
  }
  randList=pbrPar(rep(adjBlockSize,times=500),K=length(ratio),ratio=ratio,group=groupNames)
  M<-genSeq(randList,1)
  randList=as.numeric(getRandList(M))
  return(randList)
}

# addISARandomizationList
# Input: 
#	isaDesignList: the complete list of all ISA designs which have already been initiated already.
#	platformSimulationParameters: list with simulation options.
#	platformDesign: list with the overarching design parameters for the platform trial.
# Output:
#	isaDataList: updated list of ISA designs with added randomization list to each ISA
addISARandomizationList=function(isaDesignList,platformSimulationParameters,platformDesign){
	nISA=length(isaDesignList)
	for(i in seq(1,nISA)){
		if(any(names(isaDesignList[[i]])=="randomizationList")){
			isaDesignList[[i]]$randomizationList=getRandomizationListISA(isaDesignList[[i]],platformSimulationParameters,
				platformDesign)
		}
	}
	return(isaDesignList)
}

# getRandomizationListISA
# Input: 
#	isaDesign: list containing information and status of ISA
#	simulationParameters: list containing parameters used to define randomization list
#		blockSize
#	platformDesign: list containing platform design settings
# Output:
# 	Vector of assignments within ISA
getRandomizationListISA=function(isaDesign,simulationParameters,platformDesign){
	if(platformDesign$isRAR==TRUE){
		allocProb=isaDesign$allocProb
	} else {
		allocProb=isaDesign$targetSize/sum(isaDesign$targetSize)
	}
	ratio=ceiling(allocProb*simulationParameters$blockSize)
	adjBlockSize=sum(ratio)
	groupNames=as.character(seq(1,length(ratio)))
	if(any(ratio==0)){
		groupNames=groupNames[which(ratio!=0)]
		ratio=ratio[which(ratio!=0)]
	}
	randList=pbrPar(rep(adjBlockSize,times=ceiling((sum(isaDesign$targetSize)-sum(isaDesign$currentSize))/adjBlockSize)*2),K=length(ratio),ratio=ratio,group=groupNames)
	M<-genSeq(randList,1)
	randList=as.numeric(getRandList(M))
    return(randList)
}

# getVectorFromList
# Input: 
#	listIn: a list
#	name:	a string, which is a name of an element in the list
#	sumElements: boolean. If true, elements of the list element are summed up
# Output:
# 	vector as long as the list, with information on the list element "name"
getVectorFromList=function(listIn,name,sumElements=TRUE){
	out=c()
	for(i in seq(1,length(listIn))){
		index=which(names(listIn[[i]])==name)
		if(sumElements){
			out[i]=sum(listIn[[i]][[index]])
		} else {
			out[i]=listIn[[i]][[index]]
		}
	}
	return(out)
}

# listToVec
# Note: Similar to rbindlist
# Input: 
#	list: list of "length 1".
# Output:
# 	Vector containing each element of the list.
listToVec=function(list){
	out=c()
	for(i in seq(1,length(list))){
		out=c(out,list[[i]])
	}
	return(out)
}

# repRows
# Input: 
#	vec: vector of interest
#	n: number of replications
# Output:
# 	Array, which replicates vec "n" times as rows below each other
repRows=function(vec,n){
	if(length(vec)==1){
		vecOut=rep(vec[1],n)
	} else {
		vecOut=t(array(vec,dim=c(length(vec),n)))
	}
	return(vecOut)
}

# getLastVisitIndexes
# Input: 
#	data: data frame with multiple rows per patient. Visits per patient ordered.
# Output:
# 	Vector with the index position of the last visit per patient
getLastVisitIndexes=function(data){
	patID=unique(data)
	patIDVec=data
	indexVec=rep(0,length(patID))
	tmp=0
	for(i in seq(1,length(patID))){
		indexes=sum(patIDVec==patID[i])
		indexVec[i]=indexes+tmp
		tmp=indexVec[i]
	}
	out=indexVec
	return(out)
}

