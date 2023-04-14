######
### Programmed by: Tobias  Mielke
### These functions are specific to TB simulations and continue build the API to the generic simulator
### Version: Draft (run)
######
				

### functions to be customized for each platform
library(invgamma)
library(survival)
library(DescrTab2)
# interface function conducting analysis

updateTargetSize=function(allocationProbability,isaDesignListTMP){
	return(isaDesignListTMP)
}

getEnrichmentVariable=function(){
	return("CD4")
}

getControlAnalysis=function(controlISA,controlData){
	return(controlISA)
}

customizedData=function(arm,visitVec,parameters,isaDesign){
### Culture conversion
	responder=rbinom(1,size=1,prob=parameters$cultureConversionProb[arm])
	if(responder){
		lambdaPar=exp(log(-log(1-parameters$eventRate[arm]))/parameters$shapeParameter[arm])/parameters$followUp[arm]
		eventTime=rweibull(1,shape=parameters$shapeParameter[arm],scale=1/lambdaPar)
	} else {
		eventTime=Inf
		lambdaPar=NULL
	}

	visitVecTMP=visitVec
	out=array(0,dim=c(length(visitVecTMP),5))
	out[,1]=visitVecTMP
	out[,2]=rep(eventTime,length(visitVecTMP))
	
	if(any(visitVec>eventTime)){
		out[,3]=rep(visitVec[min(which(visitVec>eventTime))],length(visitVecTMP))
		out[,4]=as.numeric(visitVecTMP>=(eventTime))
	} else {
		out[,3]=Inf
		out[,4]=0
	} 
	
	unfavorableEvent=rbinom(1,size=1,prob=parameters$unfavorableEventProb[arm])
	if(eventTime>26){
		unfavorableEvent=1
	} else {
		unfavorableEvent=rbinom(1,size=1,prob=parameters$unfavorableEventProb[arm])
	}
	out[,5]=unfavorableEvent
	
	colnames(out)=c("visitID","TTCC","TTCC.1","TTCC.Censor","event")
	return(out)
	### out1 = real culture conversion time
	### out2 = detected culture conversion (visit)
	### out3 = declared culture conversion (visit+7)
	### out4 = censor

### Unfavorable event (after culture conversion)
}

customTrigger=function(isaDesign,time,data){
	nStages=length(isaDesign$interimCutoff)
	analysisCondition=rep(0,nStages)
	k=length(isaDesign$interimCutoff)
	
	if(isaDesign$activeISA){
		# Trigger for IA: Number of events on concurrent control
		subData=data[which(data$arm==0),]
		
		# Shift time by 7 weeks for culture negativity
		
		eventTime=subData$TTCC.1+7
		# Events: Prior to drop-out and event read-out prior to analysis time
		tmp1=(eventTime<subData$dropOutTime)
		tmp2=(eventTime+subData$entryTime<=time)
		nEvents=sum(tmp1*tmp2)
	
		for(i in seq(1,(k-1))){
			if(nEvents>=isaDesign$interimCutoff[i] & isaDesign$interimPassed[i]==0){
				analysisCondition[i]=i
			}
		}
	} else {
		# Final analysis triggered, when all enrolled patients reach 52 week follow-up
		# Calculated as LPI+52
		# This case is active only, when the ISA is not actively enrolling anymore.
	
		followUp=isaDesign$interimFollowUp[k]
	
		if(time>=(max(data[which(data$arm!=0),]$entryTime)+followUp)& isaDesign$interimPassed[k]==0){
			analysisCondition[k]=1
		}
	}
		if(any(analysisCondition!=0)){
			return(min(which(analysisCondition!=0)))
		} else {
			return(0)
		}
	
}

customAnalysisFunctions=function(time,data,interim,isaDesign){
	if(isaDesign$interimTrigger!="perpetual"){
		if(interim<isaDesign$firstInterim){
			tmpOut=list(coefficients=cbind(rep(0,4),rep(1,4)),varMat=diag(1,4),postProb=rep(0.6,4),nEventsCtrl=0,status=0)
			return(tmpOut)
		} else if(interim<length(isaDesign$interimCutoff)){
			tmpOut=list(coefficients=cbind(rep(0,4),rep(1,4)),varMat=diag(1,4),postProb=rep(0.6,4),nEventsCtrl=0,status=0)
			out=tryCatch(timeToCultureConversion(time,data,interim,isaDesign$calcPostProb),
				error=function(e){return(tmpOut)})
			return(out)
				
		} else {
			return(cultureConversionRate(time,data,interim))
		}
	} else {
		if(isaDesign$activeISA){
			tmpOut=list(coefficients=cbind(rep(0,4),rep(1,4)),varMat=diag(1,4),postProb=rep(0.6,4),nEventsCtrl=0,status=0)
				out=tryCatch(timeToCultureConversion(time,data,interim,isaDesign$calcPostProb),error=function(e){return(tmpOut)})
		} else {
			return(cultureConversionRate(time,data,interim))
		}
	}
}

# interface function evaluating decision
getDecisionCustom=function(interim,isaDesign,analysisResults){
	if(isaDesign$interimTrigger!="perpetual"){
		nAnalyses=length(isaDesign$interimCutoff)
		if(interim<nAnalyses){
			return(regularInterimAnalysis(isaDesign,analysisResults))
		} else {
			return(finalAnalysis(isaDesign,analysisResults))
		}
	} else {
		if(isaDesign$activeISA){	
			return(regularInterimAnalysis(isaDesign,analysisResults))
		} else {
			return(finalAnalysis(isaDesign,analysisResults))
		}
	
	}
}

# checking whether lower one-sided CI > Minimum acceptable value (MAV) (GO)
## Testing for favorable event rate => 5% vs. 15%. Success, if delta significantly below 10%
testMAV=function(data,mav,level){
	nArms=dim(data)[1]
	out=c()
	
	for(i in seq(1,nArms-1)){
		nActive=data[i+1,3]
		sActive=data[i+1,1]*nActive
		resA=c(rep(1,sActive),rep(0,nActive-sActive))
		nCtrl=data[1,3]
		sCtrl=data[1,1]*nCtrl
		resC=c(rep(1,sCtrl),rep(0,nCtrl-sCtrl))
	
		out[i]=tryCatch(farrington.manning(resA==1,resC==1,delta=mav,alternative="less",alpha=level)$p.value<level,
				error=function(e){NA})
		if(is.na(out[i])){
			return(data[i+1,1]-data[1,1]+qnorm(1-level)*sqrt(data[i+1,2]^2/data[i+1,3]+data[1,2]^2/data[1,3])<mav)}		
	}
	return(out)
}


splitData=function(assessment,data){
	assessment=sort(c(1e-10,assessment))
	timeStart=c()
	timeEnd=c()
	event=c()
	for(i in seq(1,dim(data)[1])){
		if(data$eventCensor[i]==0){
			timeStart[i]=assessment[max(which(assessment<=data$followUp[i]))]
			timeEnd[i]=Inf
			event[i]=0
		} else {
			timeStart[i]=assessment[max(which(assessment<data$followUp[i]))]
			timeEnd[i]=assessment[min(which(assessment>=data$followUp[i]))]
			if(timeStart[i]==0){
				event[i]=2
			} else {
				event[i]=3
			}
		}
	}
	outDat=data
	outDat$tstart=timeStart
	outDat$tstop=timeEnd
	outDat$eventIC=event
	return(outDat)	
}

getCorrectedFollowUp=function(data,visitVector){
	for(i in seq(1,dim(data)[1])){
		if(any(visitVector<data$TTCC.1[i])){
			lastPositiveVisit=visitVector[max(which(visitVector<data$TTCC.1[i]))]
			data$followUp[i]=max(1e-10,lastPositiveVisit)
		} else {
			data$followUp[i]=1e-10
		}
	}
	return(data)
}

getSurvivalData=function(time,data){
	visitVector=unique(data$visitID)
	indexes=getLastVisitIndexes(data$patID_2)
	data=data[indexes,]
	data$followUp=data$visitTime-data$entryTime
	data=data[which(data$followUp>0),]
	data$eventCensor=(data$entryTime+data$TTCC.1+7<=time)
	if(any(data$eventCensor==1)){
		data$followUp[which(data$eventCensor==1)]=data$TTCC.1[which(data$eventCensor==1)]
	}
	if(any((data$entryTime+data$TTCC.1)<=time & (data$entryTime+data$TTCC.1+7)>time)){
		problemCase=which((data$entryTime+data$TTCC.1)<=time & (data$entryTime+data$TTCC.1+7)>time)
		data[problemCase,]=getCorrectedFollowUp(data[problemCase,],visitVector)
	}
	
	tmpDropOut=(data$dropOutTime<=data$followUp)
	if(any(tmpDropOut)){
		data$followUp[which(tmpDropOut)]=data$dropOutTime[which(tmpDropOut)]
		data$TTCC.Censor[which(tmpDropOut)]=0
		data$eventCensor[which(tmpDropOut)]=0
	}
	data=splitData(visitVector,data)
	return(data)	
}

censorAtTrtCompletion=function(data){
	censorAtCompletion=data$trtSpecifics<data$TTCC.1
	indexes=which(censorAtCompletion)
	data$tstart[indexes]=data$trtSpecifics[indexes]
	data$tstop[indexes]=Inf
	data$eventIC[indexes]=0
	return(data)
}

simulateTime=function(data,shape,scale){
	data$simTime=data$tstart
	if(sum(data$eventIC==3)>0){
		indexes=which(data$eventIC==3)
		tStart=data$tstart[indexes]
		tStop=data$tstop[indexes]
		Y=runif(sum(data$eventIC==3),min=0,max=1)
		tmp=exp(-(tStart/scale)^shape)-Y*(exp(-(tStart/scale)^shape)-exp(-(tStop/scale)^shape))
		out=scale*(-log(tmp))^(1/shape)
		data$simTime[indexes]=out
	}
return(data)
}


timeToCultureConversion=function(time,data,interim,calcPostProb){
	data=getSurvivalData(time,data)
	nArms=length(unique(data$arm))
	
	### Weibull Analysis
	fullAna=summary(survreg(Surv(time=data$tstart,
		time2=data$tstop,
		event=eventIC,
		type="interval")~1+as.factor(arm),dist="weibull",data=data))

	coefOut=cbind(fullAna$coefficients[-c(1,nArms+1)],sqrt(diag(fullAna$var)[-c(1,nArms+1)]))
	varOut=fullAna$var[-(nArms+1),-(nArms+1)]
	varOut=varOut[-1,-1]
	
	ctrlDat=data[data$arm==0,]
	nEventsCtrl=sum(ctrlDat$eventIC==3)
	
	if(calcPostProb){
		tmpCtrl=survreg(Surv(time=ctrlDat$tstart,
			time2=ctrlDat$tstop,
			event=eventIC,
			type="interval")~1,dist="weibull",data=ctrlDat)

		shapeCtrl=1/exp(tmpCtrl$icoef[2])

		xCtrl=c()
		for(i in seq(1,10)){
			datTMP=simulateTime(ctrlDat,shapeCtrl,scale=exp(coef(tmpCtrl)))
			aCtrl=sum(datTMP$eventIC==3)
			bCtrl=sum(datTMP$simTime^shapeCtrl)
			xCtrl=c(xCtrl,rinvgamma(1000,shape=aCtrl,rate=bCtrl)^(1/shapeCtrl))
		}
		
		activeDat=data[data$arm!=0,]
		activeDat=censorAtTrtCompletion(activeDat)
		tmpAct=survreg(Surv(time=activeDat$tstart,
			time2=activeDat$tstop,
			event=eventIC,
			type="interval")~1,dist="weibull",data=activeDat)
		shapeActive=1/exp(tmpAct$icoef[2])

		xActive=c()
		for(i in seq(1,10)){
			datTMP=simulateTime(ctrlDat,shapeActive,scale=exp(coef(tmpAct)))
			aActive=sum(datTMP$eventIC==3)
			bActive=sum(datTMP$simTime^shapeActive)
			xActive=c(xActive,rinvgamma(1000,shape=aActive,rate=bActive)^(1/shapeActive))
		}
		
		trtSpecifics=sort(unique(data$trtSpecifics))
		trtSpecifics=c(12,16,20,26)
		probCtrl=1-exp(-(26/xCtrl)^shapeCtrl)
		
		probActive=c()
		postProb=c()
		
		niDelta=0.1
		mnAct=c()
		for(i in seq(1,length(trtSpecifics))){
			probActive=1-exp(-(trtSpecifics[i]/xActive)^shapeActive)
			mnAct[i]=mean(probActive)
			postProb[i]=mean(probActive>(probCtrl-niDelta))
		}
	} else {
		trtSpecifics=sort(unique(data$trtSpecifics))
		trtSpecifics=c(12,16,20,26)
		postProb=rep(1,length(trtSpecifics))
	}
	out1=list(coefficients=coefOut,varMat=varOut,postProb=postProb,nEventsCtrl=nEventsCtrl,status=1)
}

# actual function getting the relevant data summary
cultureConversionRate=function(time,data,interim){
		subData=data

		subData=data[visitID==52,]
		nArms=length(unique(subData$arm))
		
		rate=c()
		n=c()
		duration=c()
		uArms=sort(unique(subData$arm))
		for(i in seq(1,nArms)){
			indexes=which(subData$arm==uArms[i])
			rate[i]=mean(subData$event[indexes])
			n[i]=length(indexes)
			duration[i]=subData$trtSpecifics[indexes[1]]
		}

		out2=as.matrix(cbind(rate,n,duration))
		return(out2)
}

# function getting some statistics for binary data (relative frequency, sample size, standard deviation)
ratesToNormal=function(data){
	out=array(0,dim=c(dim(data)[1],3))
	out[,1]=data[,1]
	out[,3]=data[,2]
	out[,2]=sqrt(data[,1]*(1-data[,1]))
	return(out)
}

### summary
getDecisionMat=function(isaDesign){
	
	if(isaDesign$decisionMaking$interimTrigger=="perpetual"|isaDesign$decisionMaking$interimTrigger=="customTrigger"){
		nInterim=1
	} else {
		nInterim=length(isaDesign$decisionMaking$interimCutoff)
	}
	nArms=length(isaDesign$designElements$targetSize)
	colNames=c("Rate_Active_","HR_","Go_arm_","Success_arm_","Size_arm_")
	colNames=paste(rep(colNames,times=nArms),rep(seq(1,nArms),each=length(colNames)),sep="")
	
	out=array(0,dim=c(nInterim,length(colNames)+5))
	colnames(out)=c("Stage","SummaryGo","Control_Rate","Size_ctrl",colNames,"Final")
	return(out)
}

getAdaptiveCutoffs=function(interim,isaDesign,sampleSize,decision,analysisResults){
	cutoffsOut=isaDesign$interimCutoff
	if(isaDesign$interimTrigger=="perpetual"){
			return(cutoffsOut)
	}
	if(isaDesign$interimTrigger=="events"){
			return(cutoffsOut)
	} else if(isaDesign$interimTrigger=="time"){
			return(cutoffsOut)
	} else {
			if(any(sampleSize<isaDesign$targetSize)){
				cutoffsOut[length(isaDesign$interimCutoff)]=sum(sampleSize)
			}
			return(cutoffsOut)
	}
}		

regularInterimAnalysis=function(isaDesign,analysisResults){
	isa=isaDesign[[2]]
	stage=isaDesign$stage
	nArms=dim(analysisResults[[1]][[1]])[1]-1
	nArms=length(isaDesign$targetSize)
	
	postProb=analysisResults[[1]][[3]]
	summary=rep(0,length=length(colnames(isaDesign$decisions)))

	summary[1]=stage
	summary[3]=NA
	
	rateColumns=seq(5,5+(nArms-1)*5,by=5)
	summary[rateColumns]=NA
	
	sizeColumns=seq(9,9+(nArms-1)*5,by=5)
	summary[sizeColumns]=isaDesign$currentSize
	summary[4]=isaDesign$controlSize
	
	hrColumns=seq(6,6+(nArms-1)*5,by=5)
	goColumns=seq(7,7+(nArms-1)*5,by=5)
	successColumns=seq(8,8+(nArms-1)*5,by=5)
	if(analysisResults[[1]]$status==0){
		# Analysis not successful => continue with all.
			if(stage>1){
				go=1*isaDesign$decisions[stage-1,1]
			} else {
				go=1
			}
			summary[2]=go
			summary[goColumns]=go
			summary[successColumns]=0
			cutoffs=isaDesign$interimCutoff
			population=1
			arms=rep(1,nArms)
			allocProb=isaDesign$allocProb
			
			return(list(summary=summary,sampleSize=isaDesign$targetSize,cutoffs=cutoffs,population=population,arms=arms,allocProb=allocProb))	
	} else {
		pValue=1-pnorm(-analysisResults[[1]][[1]][,1]/analysisResults[[1]][[1]][,2])
		hr=exp(-analysisResults[[1]][[1]][,1])
		
		if(any(is.na(pValue))){
			print(analysisResults[[1]])
		}
		
		summary[hrColumns]=hr
		
		informationFraction=analysisResults[[1]][[4]]/124
		
		futilityBoundary=futilityRule(informationFraction=informationFraction,
										isaDesign$decisionRules$cutoffsIF,
										isaDesign$decisionRules$cutoffsAlpha)
		
		
		if(min(pValue[which(!is.na(pValue))])>futilityBoundary){
			# All arms futile based on model-free approach
			summary[2]=0
			summary[goColumns]=0
			summary[successColumns]=0
						
			sampleSize=isaDesign$currentSize
			cutoffs=getAdaptiveCutoffs(interim=stage,isaDesign,sampleSize,decision=0,analysisResults)
			arms=rep(0,nArms)
			population=0
			return(list(summary=summary,sampleSize=sampleSize,cutoffs=cutoffs,population=population,arms=arms,allocProb=postProb))
		} else {
			# Some arms non-futile
			sampleSize=isaDesign$targetSize
			arms=rep(1,nArms)
			index=1
			
			### Indicator vector
			nonFutile=rep(0,nArms)
			nonFutileMB=rep(1,nArms)
			
			### Standard case: model free analysis
			if(length(pValue)>1){
					# What is the shortest duration, which is not stopped for futility?
					index=min(which(pValue<futilityBoundary))
			} else {
					index=1
			}
			nonFutile[index:nArms]=1
			
			if(any(isaDesign$decisionRules$postProbCutoff>1e-5)){			
				# Do this only, in case that posterior probability is used
				# Otherwise, nonFutileMB is 1 for all arms.
				
				armFutilityBoundary=futilityRule(informationFraction=informationFraction,
										isaDesign$decisionRules$cutoffsIF,
										isaDesign$decisionRules$postProbCutoff)
				if(any(postProb<armFutilityBoundary)){
				 indexes=which(postProb<armFutilityBoundary)
				 nonFutileMB[indexes]=0
				}
			}
		
			arms=nonFutile*nonFutileMB
			indexes=which(arms==0)
			sampleSize[indexes]=isaDesign$currentSize[indexes]
					
			if(stage>1){
				go=1*isaDesign$decisions[stage-1,1]
			} else {
				go=1
			}
			summary[2]=go
			summary[goColumns]=arms
			summary[successColumns]=0
			if(sum(goColumns)==0){
				summary[2]=0
			}

			cutoffs=getAdaptiveCutoffs(interim=stage,isaDesign,sampleSize,decision=1,analysisResults)
			population=1
			
			return(list(summary=summary,sampleSize=sampleSize,cutoffs=cutoffs,population=population,arms=arms,allocProb=postProb))	
		}
	}
}

finalAnalysis=function(isaDesign,analysisResults){

	nArms=length(isaDesign$targetSize)
	stage=isaDesign$stage
	analysisResultsTMP=analysisResults[[1]][1:2,]	
	
	summary=rep(0,length=length(colnames(isaDesign$decisions)))
	summary[1]=stage
	summary[2]=NA
	
	rateColumns=seq(5,5+(nArms-1)*5,by=5)
	summary[rateColumns]=analysisResults[[1]][-1,1]
	summary[3]=analysisResults[[1]][1,1]
	sizeColumns=seq(9,9+(nArms-1)*5,by=5)
	summary[sizeColumns]=isaDesign$currentSize
	summary[4]=isaDesign$controlSize

	hrColumns=seq(6,6+(nArms-1)*5,by=5)
	summary[hrColumns]=NA
		
	goColumns=seq(7,7+(nArms-1)*5,by=5)
	successColumns=seq(8,8+(nArms-1)*5,by=5)
	
	successComplete=rep(0,nArms)
	successFinal=rep(0,nArms)
	finalMAV=isaDesign$decisionRules$MAV[length(isaDesign$decisionRules$MAV)]

	for(i in seq(1,nArms)){
		# columns: response rate, samples size, trt duration.
		if(any(analysisResults[[1]][,3]==isaDesign$trtSpecifics[i])){
			# This is case is to cover the (unlikely situation) that no patient has been randomized to one of the arms.
			arm=which(analysisResults[[1]][,3]==isaDesign$trtSpecifics[i])
			if(isaDesign$trtSpecifics[i]==26){
				arm=arm[2]
			}
	
			analysisResultsTMP[2,]=c(sum(analysisResults[[1]][arm,1]*analysisResults[[1]][arm,2]),sum(analysisResults[[1]][arm,2]),24)
			analysisResultsTMP[2,1]=analysisResultsTMP[2,1]/analysisResultsTMP[2,2]
			mav=testMAV(ratesToNormal(analysisResultsTMP),finalMAV,isaDesign$decisionRules$levelFinal)
			## Success for the arm: Final analysis success and not stopped in IA
			successComplete[i]=(mav*prod(isaDesign$decisions[(1:(stage-1)),goColumns[i]]))>0
			successFinal[i]=mav
		} else {
			successComplete[i]=0
			successFinal[i]=0
		} 
	}
	summary[1]=999
	summary[2]=successComplete[nArms]
	summary[successColumns]=successComplete
	summary[goColumns]=successFinal
	summary[length(summary)]=1
	return(list(summary=summary,sampleSize=isaDesign$currentSize,cutoffs=NA,population=0,arms=rep(0,nArms),allocProb=rep(0,nArms)))
}


getEffectScenario=function(relativeEffect,worstCase,bestCase,shapeParameter,followUp){
	out=array(0,dim=c(length(relativeEffect),length(worstCase)))

	for(i in seq(1,length(relativeEffect))){
		out[i,]=worstCase+relativeEffect[i]*(bestCase-worstCase)
	}
	
	return(list(
		cultureConversionProb=out[,1],
		eventRate=out[,2],
		followUp=rep(followUp,times=length(relativeEffect)),
		shapeParameter=rep(shapeParameter,times=length(relativeEffect)),
		unfavorableEventProb=out[,3]
	))
}


getCustomizedAllocationProbability=function(isaDesignList,platformDesign){
	nISA=length(isaDesignList)
	for(i in seq(1,nISA)){
		openEnrollmentPerISA=getVectorFromList(isaDesignList,"activeISA",sumElements=TRUE)
	}
	nActiveArms=c(0,0)
	allocProb=rep(0,nISA)
	## nNewISA: Vector indicating, that ISA had no interim analysis. Otherwise, take "run-in" allocation.
	nNewISA=rep(0,nISA)
	## oldOpenISA: Vector indicating, whether ISA already had interim analysis and is active
	## mean of allocation accross probability across ISAs, which had IA will be used for new ISA to start.
	oldOpenISA=rep(0,nISA)
	tmpAlloc=c()
	for(i in seq(3,nISA)){
		if(openEnrollmentPerISA[i]){
			activeArms=isaDesignList[[i]]$currentSize<isaDesignList[[i]]$targetSize
			nActiveArms[i]=sum(activeArms)
			if(isaDesignList[[i]]$allocationRule=="constant"){
				allocProb[i]=sqrt(isaDesignList[[i]]$allocProb[2])
				if(allocProb[i]<0.05){
					allocProb[i]=0.05
				}
				tmpAlloc[i]=isaDesignList[[i]]$allocProb[2]
			} else {
				allocProb[i]=mean(sqrt(isaDesignList[[i]]$allocProb))
				tmpAlloc[i]=mean(isaDesignList[[i]]$allocProb)
				if(allocProb[i]<0.05){
					allocProb[i]=0.05
				}
			}
			nNewISA[i]=sum(isaDesignList[[i]]$analysisTime>0)==1
			oldOpenISA[i]=sum(isaDesignList[[i]]$analysisTime>0)>1
		} else {
			nActiveArms[i]=0
			allocProb[i]=0
			nNewISA[i]=0
			oldOpenISA[i]=0
		}
	}
	if(sum(nNewISA>0)){
		if(sum(oldOpenISA)>0){
				allocProbTMP=sqrt(mean(tmpAlloc[which(oldOpenISA==1)]))
		} else {
				allocProbTMP=1
		}
		allocProb[which(nNewISA==1)]=allocProbTMP
	}
	# For RAR in this simulation only use square-root allocation to control
	ctrlAlloc=1/(sqrt(sum(nActiveArms))+1)
	ctrlBoundary=as.numeric(platformDesign$allocationRule[3])
	if(ctrlAlloc<ctrlBoundary){
		ctrlAlloc=ctrlBoundary
	}
		
	allocProbTMP=allocProb
	if(sum(allocProb)==0){
		allocProb[2]=1
	} else {
		allocProb=allocProb/sum(allocProb)
		allocProb=c(0,ctrlAlloc,(1-ctrlAlloc)*allocProb[-c(1:2)])
	}
	return(allocProb)	
}

