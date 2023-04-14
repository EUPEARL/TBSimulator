######
### Programmed by: Tobias  Mielke
### These functions are specific to TB for the definition of the ISA designs
### Version: Draft (run)
######
									

futilityRule=function(informationFraction,cutoffsIF,cutoffsAlpha){
	if(informationFraction<=cutoffsIF[1]){
	 out=cutoffsAlpha[1]
	} else {
		if(any(informationFraction<=cutoffsIF)){
			index=min(which(informationFraction<=cutoffsIF))
			if(index<=length(cutoffsIF)){
			out=cutoffsAlpha[index-1]+(informationFraction-cutoffsIF[index-1])/
				(cutoffsIF[index]-cutoffsIF[index-1])*(cutoffsAlpha[index]-cutoffsAlpha[index-1])
			}
		} else {
			out=cutoffsAlpha[length(cutoffsAlpha)]
		}
	}
	return(out)
}


getDesignElements=function(nArms,trtDuration,priorAllocProb,sizePerArm,dropOutRate,allocationRule){
	designElements=list(
		targetSize=rep(sizePerArm,nArms),
		plannedSize=rep(sizePerArm,nArms),
		minSize=rep(sizePerArm,nArms),
		maxSize=rep(sizePerArm,nArms),
		visitVec=c(0,2,4,6,8,12,17,26,39,52),
		randomizationRatio=rep(1,nArms),
		allocProb=priorAllocProb,
		trtSpecifics=trtDuration,
		allocationRule=allocationRule,
		dropOutEventRate=rep(dropOutRate,nArms))
	return(designElements)
}

getDecisionRules=function(ia,postProbCutoff,sizePerArm){
calcPostProb=min(postProbCutoff>=0)

levelFE=0.1
niMargin=0.1

if(ia==4){
	decisionRules_4IA=list(
		MAV=c(1,1,1,1,niMargin),
		cutoffsAlpha=c(0.8,0.4,0.2,0.2),
		cutoffsIF=c(10/124,28/124,50/124,75/124),
		levelFinal=levelFE,
		postProbCutoff=postProbCutoff)
	
	decisionMaking=list(
		firstInterim=1,
		interimTrigger="customTrigger",
		interimCutoff=c(10,28,50,75,sizePerArm*4),
		interimFollowUp=c(12,12,12,12,52),
		decisionRules=decisionRules_4IA,
		calcPostProb=calcPostProb)
} else if(ia==2){
	decisionRules_2IA=list(
		MAV=c(1,1,1,1,niMargin),
		cutoffsAlpha=c(1,1,0.4,0.2),
		cutoffsIF=c(5/124,10/124,28/124,50/124),
		
		levelFinal=levelFE,
		postProbCutoff=postProbCutoff)

	decisionMaking=list(
		firstInterim=3,
		interimTrigger="customTrigger",
		interimCutoff=c(5,6,28,50,sizePerArm*4),
		interimFollowUp=c(12,12,12,12,52),
		decisionRules=decisionRules_2IA,
		calcPostProb=calcPostProb
		)
}
	return(decisionMaking)
}

getAnalysisElements=function(){
analysisElements=list(
	nEndpoints=2,
	controlType="concurrentControl",
	modelBased=FALSE)
return(analysisElements)
}

getDashboard=function(ia,activeISA){
	if(ia==2){
		nRep=5
	} else if(ia==4){
		nRep=5
	} else {
		nRep=2
	}
	dashboard=list(
		stage=1,
		activeISA=activeISA,
		completed=0,
		currentSize=0,
		controlSize=0,
		interimPassed=rep(0,nRep),
	decisions=NA)
}

getDesign=function(nArms,sizePerArm,trtDuration,priorAllocProb,postProbCutoff,ia,activeISA,name,scenario,dropOutRate,allocationRule){
	designElements=getDesignElements(nArms,trtDuration,priorAllocProb,sizePerArm,dropOutRate,allocationRule)
	decisionRules=getDecisionRules(ia,postProbCutoff,sizePerArm)
	analysisElements=getAnalysisElements()
	dashboard=getDashboard(ia,activeISA)
	
	outList=list(
		name=name,
		effectScenario=scenario,
		designElements=designElements,
		analysisElements=analysisElements,
		decisionMaking=decisionRules,
		dashboard=dashboard)
	return(outList)
}

defineISAList=function(isaList,isaIndex){
	isaListOut=list()
for(i in seq(1,length(isaIndex))){
	isaListOut[[i]]=isaList[[isaIndex[i]]]
	isaListOut[[i]]$name=paste("ISA_",i-2,sep="")
	isaListOut[[i]]$effectScenario=i
}
isaListOut[[1]]$name="long"
isaListOut[[1]]$effectScenario=1
isaListOut[[2]]$name="control"
isaListOut[[2]]$effectScenario=2
return(isaListOut)
}

### Rec-separate is step function:
# Element 1: Recruitment rate at t=0
# Element 2: Time point with change in rate
# Element 3: Additional rec rate from time-point 2 onwards
# Element 4: Timepoint from which only element 3 used as rate.

### Simulation settings
# First element: Number of inactive weeks, when platform is terminated
# Second element: Interval of monitoring for IA
# Third element: Maximum data size prior to merge with full data

getSimulationParameters=function(recruitmentSettings,baselineSimulationParameters,computationalSettings,endpointSettings){
simPars=list()
for(i in seq(1,length(recruitmentSettings))){
	simPars[[i]]=list(
			recruitmentPars=recruitmentSettings[[i]]$recruitmentParameters,
			recruitmentType=recruitmentSettings[[i]]$recruitmentType,
			baselineSimulationPars=baselineSimulationParameters,
			randomizationParameter=computationalSettings[[1]],
			blockSize=computationalSettings[[2]],
			useRandomizationList=computationalSettings[[3]],
			simulationSettings=computationalSettings[[4]],
			storeAllocRates=computationalSettings[[5]],
			nEndpoints=endpointSettings[[1]],
			endpoint=endpointSettings[[2]])
	}
return(simPars)
}


getPlatformDesign=function(allocProb,allocRule,entryTime,
					isaOrder,isaProb,randomizer,isaList){
	platformDesign=list(allocationRule=allocRule,
				randomizer=randomizer,
				entryTime=entryTime,
				fixedISAOrder=isaOrder,
				lateISAList=isaList,
				isaProb=isaProb)
	return(platformDesign)
}