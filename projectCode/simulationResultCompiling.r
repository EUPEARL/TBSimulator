######
### Programmed by: Tobias  Mielke
### These functions are specific to some TB simulations to post-process output from the platform trial simulations
### Version: Draft (run)
######
				

fdr=function(rejections,type1){
	if(sum(rejections)==0){
		return(0)
	} else {
		return(sum(type1)/sum(rejections))
	}
}

getNRejections=function(out,nISA){
mat=cbind(out$successArm1,out$successArm2,out$successArm3,out$successArm4)
sum=0
for(i in seq(1,nISA-2)){
#for(i in seq(1,2)){
	if(mat[i,4]==1){
		if(mat[i,3]==1){
			if(mat[i,2]==1){
				if(mat[i,1]==1){
					sum=sum+4
				} else {
					sum=sum+3
				}
			} else {
					sum=sum+2
			}
		} else {
			sum=sum+1
		}
	}
}
return(sum)					
}
# Success = at least one intersection (corrected)
# SuccesArm1 = number intersection(corrected)
# SuccessArm2 = number total (corrected)
# SuccessArm3 = number intersection (non corrected)
# SuccessArm4 = number total (non corrected)

## need false discover rate on arms
## number rejected corrected
## number type-1 errors corrected.
## Ordering of results an issue

getCtrlSummary=function(ctrlISA,out,nISA){
	ctrlISA=c(2,1,max(out$Stop),max(out$Stop),rep(NA,8),max(out$atLeastOneSuccess),sum(out$atLeastOneSuccess),rep(NA,3),sum(out$nSuccess),
	max(out$atLeastOneSuccessNC),sum(out$atLeastOneSuccessNC),rep(NA,3),sum(out$nSuccessNC),rep(NA,5),
	sum(out$totalActiveSize)+ctrlISA$currentSize,sum(out$totalActiveSize),ctrlISA$currentSize,rep(NA,4),
	max(out$atLeastOneType1),sum(out$atLeastOneType1),rep(NA,2),ifelse(sum(out$nSuccess)==0,0,sum(out$nType1)/sum(out$nSuccess)),sum(out$nType1),
	max(out$atLeastOneNC),sum(out$atLeastOneNC),rep(NA,2),ifelse(sum(out$nSuccessNC)==0,0,sum(out$nType1NC)/sum(out$nSuccessNC)),sum(out$nType1NC),
	mean(out$rate0),rep(NA,4))
}

### Success statistics
# max(success) = at least one (corrected for mult arm 4)
# sum(success) = total number of successes (corrected for mult - arm 4)
# getNRejections = total number of success (corrected for all arms)
# sum(out$nonCorrected) = at least one (if not corrected for multiplicity)
# sum(out$successArm1+...) = total number of successes (not corrected for mult)
###

### Type-1 statistics
# max(out$Type1_Corrected)= at least one type-1 error (arm 4)
# fdr(out$success,out$Type1_Corrected) = at least one type-1 error (arm 4)
# fdr(out$nonCorrected,out$Type1_nonCorrected) = at least one type-1 error (non corrected = maximum over the 4 arms) - Not as targeted!
# max(out$Type1_nonCorrected) = at least one type-1 error (if you don't correct)

# Missing 1: Type-1 error if you follow the sequence of testing
# Missing 2: FDR as number of falsely rejected hypotheses among rejected hypotheses, if you were to correct (ALL)
# Missing 3: FDR as number of falsely rejected hypotheses among rejected hypotheses, if you weren't to correct (ALL)


### ctrlSummary
# ISA
# Start
# Stop
# Duration
# Stop per arm, duration per arm,  8
# At least one success
# Number of successes (Corrected)
# Total number of success (corrected)
# Number of success non corrected
# Total number of success (non corrected)
# 6 empty 
# Total size
# Total active size
# Ctrl size
# 4 empty
# FWER corrected
# FDR corrected
# FDR non-corrected
# NA
# FWER (non corrected)


getSummaryISA=function(isa,out,scenarioSummary){
	return(c(getDurationISA(isa,out),
			getSuccessISA(isa,out),
			getFutilityISA(isa,out),
			getSizeISA(isa,out),
			getType1ErrorISA(isa,out,scenarioSummary),
			getResponseRateISA(isa,out)))
}

getType1ErrorISA=function(isa,out,scenarioSummary){
success=out[(Stage==999),.(
	atLeastOneType1=(Success_arm_4*(scenarioSummary[4]+Success_arm_3*(scenarioSummary[3]+Success_arm_2*(scenarioSummary[2]+Success_arm_1*scenarioSummary[1])))>0),
	type1Arm1=Success_arm_4*Success_arm_3*Success_arm_2*Success_arm_1*scenarioSummary[1],
	type1Arm2=Success_arm_4*Success_arm_3*Success_arm_2*scenarioSummary[2],
	type1Arm3=Success_arm_4*Success_arm_3*scenarioSummary[3],
	type1Arm4=Success_arm_4*scenarioSummary[4],
	nType1=Success_arm_4*(scenarioSummary[4]+Success_arm_3*(scenarioSummary[3]+Success_arm_2*(scenarioSummary[2]+Success_arm_1*scenarioSummary[1]))),
	atLeastOneNC=max(Success_arm_1*scenarioSummary[1],Success_arm_2*scenarioSummary[2],Success_arm_3*scenarioSummary[3],Success_arm_4*scenarioSummary[4]),
	type1Arm1NC=Success_arm_1*scenarioSummary[1],
	type1Arm2NC=Success_arm_2*scenarioSummary[2],
	type1Arm3NC=Success_arm_3*scenarioSummary[3],
	type1Arm4NC=Success_arm_4*scenarioSummary[4],
	nType1NC=Success_arm_1*scenarioSummary[1]+Success_arm_2*scenarioSummary[2]+Success_arm_3*scenarioSummary[3]+Success_arm_4*scenarioSummary[4])]
out=as.numeric(success)
names(out)=c("atLeastOneType1","Type1_Arm1","Type1_Arm2","Type1_Arm3","Type1_Arm4","nType1","atLeastOneNC","Type1_Arm1_NC","Type1_Arm2_NC","Type1_Arm3_NC","Type1_Arm4_NC","nType1NC")
return(out)
}


getSizeISA=function(isa,out){
size=out[(Stage==999),.(totalSize=Size_ctrl+Size_arm_1+Size_arm_2+Size_arm_3+Size_arm_4,
				totalActiveSize=Size_arm_1+Size_arm_2+Size_arm_3+Size_arm_4,
				crltSize=Size_ctrl,
				sizeArm1=Size_arm_1,
				sizeArm2=Size_arm_2,
				sizeArm3=Size_arm_3,
				sizeArm4=Size_arm_4)]
out=as.numeric(size)
names(out)=c("totalSize","totalActiveSize","ctrlSize","sizeArm1","sizeArm2","sizeArm3","sizeArm4")
return(out)
}

### Change to include corrected success
getSuccessISA=function(isa,out){
success=out[(Stage==999),.(
	atLeastOneSuccess=(Success_arm_4*(1+Success_arm_3*(1+Success_arm_2*(1+Success_arm_1*1)))>0),
	successArm1=Success_arm_1*Success_arm_2*Success_arm_3*Success_arm_4,
	successArm2=Success_arm_2*Success_arm_3*Success_arm_4,
	successArm3=Success_arm_3*Success_arm_4,
	successArm4=Success_arm_4,
	nSuccess=(Success_arm_4*(1+Success_arm_3*(1+Success_arm_2*(1+Success_arm_1*1)))),
	atLeastOneSuccessNC=max(Success_arm_1,Success_arm_2,Success_arm_3,Success_arm_4),
	successArm1NC=Success_arm_1,
	successArm2NC=Success_arm_2,
	successArm3NC=Success_arm_3,
	successArm4NC=Success_arm_4,
	nSuccessNC=(Success_arm_4+Success_arm_3+Success_arm_2+Success_arm_1))]
	
out=as.numeric(success)
names(out)=c("atLeastOneSuccess","successArm1","successArm2","successArm3","successArm4","nSuccess","atLeastOneSuccessNC","successArm1_NC","successArm2_NC","successArm3_NC","successArm4_NC","nSuccessNC")
return(out)
}

getFutilityISA=function(isa,out){
nStage=dim(out)[1]-1
earlyStop=out[(Stage!=999),.(fut_1=sum(Go_arm_1)<nStage,
	fut_2=sum(Go_arm_2)<nStage,
	fut_3=sum(Go_arm_3)<nStage,
	fut_4=sum(Go_arm_4)<nStage,
	fut_All=(sum(Go_arm_1)<nStage)*(sum(Go_arm_2)<nStage)*(sum(Go_arm_3)<nStage)*(sum(Go_arm_4)<nStage))]
out=c(as.numeric(earlyStop))
#,nStage*as.numeric(earlyStop)[5])
names(out)=c("fut_1","fut_2","fut_3","fut_4","fut_All")
#,"stage")
return(out)
}
getAnaIndex=function(x){
	if(any(x==0)){
		return(min(which(x==0)))
	} else {
		return(length(x))
	}
}

getResponseRateISA=function(isa,out){
finalAna=as.numeric(out[(Stage==999),.(rate_0=Control_Rate,
	rate_1=Rate_Active_1,
	rate_2=Rate_Active_2,
	rate_3=Rate_Active_3,
	rate_4=Rate_Active_4)])

out=as.numeric(finalAna)
names(out)=c("rate0","rate1","rate2","rate3","rate4")
return(out)
}

getDurationISA=function(isa,out){
finalAna=as.numeric(out[,.(ana_Arm_1=getAnaIndex(Go_arm_1),
	ana_Arm_2=getAnaIndex(Go_arm_2),
	ana_Arm_3=getAnaIndex(Go_arm_3),
	ana_Arm_4=getAnaIndex(Go_arm_4))])

anaTime=out$tAna[finalAna]
durationTime=out$tAna[finalAna]-out$tStart[finalAna]
anaISA=max(anaTime)
durationISA=max(durationTime)

out=as.numeric(c(isa,out$tStart[1],anaISA,durationISA,anaTime,durationTime))
names(out)=c("ISA","Start","Stop","Duration","Stop_1","Stop_2","Stop_3","Stop_4","Duration_1","Duration_2","Duration_3","Duration_4")
return(out)
}

summarizeSingleSimulation=function(run,outList,data,allocRates,scenarioSummary){
	nISA=length(outList)
	anaTime=outList[[3]]$analysisTime[-1]
	summaryTMP=data.table(cbind(3,tStart=outList[[3]]$analysisTime[1],tAna=anaTime[which(anaTime!=0)],outList[[3]]$decisions[which(outList[[3]]$decisions[,1]!=0),]))
	out=getSummaryISA(3,summaryTMP,scenarioSummary[3,])
	if(nISA>3){
		for(i in seq(4,nISA)){
			anaTime=outList[[i]]$analysisTime[-1]
			summaryTMP=data.table(cbind(i,tStart=outList[[i]]$analysisTime[1],tAna=anaTime[which(anaTime!=0)],outList[[i]]$decisions[which(outList[[i]]$decisions[,1]!=0),]))
			out=rbind(out,getSummaryISA(i,summaryTMP,scenarioSummary[i,]))
		}
	}
	out=rbind(out,getCtrlSummary(outList[[2]],data.table(out),nISA))
	out=cbind(run=run,out)
return(out)
}

getDesignCharacteristics=function(results,path,scenario,isaList){
if(is.null(results)){
	out=read.table(paste(path,scenario,".csv",sep=""),sep="\t",header=TRUE)
	out=data.table(out)
} else {
	out=results
}

	nISA=max(unique(out$ISA))
	nSim=dim(out)[1]/(nISA-1)

	designNames=c()
	for(i in seq(3,nISA)){
		designNames[i-2]=isaList[[i]]$name
	}
	designNames[nISA-1]=isaList[[2]]$name
	out=cbind(out,design=rep(designNames,nSim))

	outputSize=out[,.(power=mean(atLeastOneSuccess),
		asn=mean(totalSize),
		nNoControl=mean(totalActiveSize),
		ctrlSize=mean(ctrlSize),
		nArm1=mean(sizeArm1),
		nArm2=mean(sizeArm2),
		nArm3=mean(sizeArm3),
		nArm4=mean(sizeArm4)),
		by=.(ISA)]
	outputPower=out[,.(power=mean(atLeastOneSuccess),
		power1=mean(successArm1),
		power2=mean(successArm2),
		power3=mean(successArm3),
		power4=mean(successArm4),
		nSuccess=mean(nSuccess),
		powerNC=mean(atLeastOneSuccessNC),
		power1NC=mean(successArm1_NC),
		power2NC=mean(successArm2_NC),
		power3NC=mean(successArm3_NC),
		power4NC=mean(successArm4_NC),
		nSuccessNC=mean(nSuccessNC)),by=.(ISA)]	
	outputType1=out[,.(fwer=mean(atLeastOneType1),
		type1Arm1=mean(Type1_Arm1),
		type1Arm2=mean(Type1_Arm2),
		type1Arm3=mean(Type1_Arm3),
		type1Arm4=mean(Type1_Arm4),
		nType1=mean(nType1),
		fwerNC=mean(atLeastOneNC),
		type1Arm1NC=mean(Type1_Arm1_NC),
		type1Arm2NC=mean(Type1_Arm2_NC),
		type1Arm3NC=mean(Type1_Arm3_NC),
		type1Arm4NC=mean(Type1_Arm4_NC),
		nType1NC=mean(nType1NC)),by=.(ISA)]
	outputFutility=out[,.(futility=mean(fut_All),
		fut1=mean(fut_1),
		fut2=mean(fut_2),
		fut3=mean(fut_3),
		fut4=mean(fut_4)),
		#stage1=mean(stage==1),
		#stage2=mean(stage==2),
		#stage3=mean(stage==3),
		#stage4=mean(stage==4)),
		by=.(ISA)]
	outputRate=out[,.(
		rateControl=mean(rate0),
		rate1=mean(rate1),
		rate2=mean(rate2),
		rate3=mean(rate3),
		rate4=mean(rate4)),by=.(ISA)]
	outputDuration=out[,.(
		duration=mean(Duration),
		start=mean(Start),
		arm1=mean(Duration_1),
		arm2=mean(Duration_2),
		arm3=mean(Duration_3),
		arm4=mean(Duration_4)),by=.(ISA)]
	return(list(sampleSize=outputSize,power=outputPower,futility=outputFutility,responseRate=outputRate,duration=outputDuration,outputType1=outputType1))
}

compareResults=function(results,isaScenarioList,designParams){
  k=1
  for(i in seq(1,length(isaScenarioList))){
    for(j in seq(1,length(designParams))){
      if(k==1){
        sampleSizeMat=cbind(scenario=i,design=j,results[[k]]$sampleSize)
        powerMat=cbind(scenario=i,design=j,results[[k]]$power)
        futilityMat=cbind(scenario=i,design=j,results[[k]]$futility)
        responseRateMat=cbind(scenario=i,design=j,results[[k]]$responseRate)
        durationMat=cbind(scenario=i,design=j,results[[k]]$duration)
		type1Mat=cbind(scenario=i,design=j,results[[k]]$outputType1)
      } else {
        sampleSizeMat=rbind(sampleSizeMat,cbind(scenario=i,design=j,results[[k]]$sampleSize))
        powerMat=rbind(powerMat,cbind(scenario=i,design=j,results[[k]]$power))
        futilityMat=rbind(futilityMat,cbind(scenario=i,design=j,results[[k]]$futility))
        responseRateMat=rbind(responseRateMat,cbind(scenario=i,design=j,results[[k]]$responseRate))
        durationMat=rbind(durationMat,cbind(scenario=i,design=j,results[[k]]$duration))
		type1Mat=rbind(type1Mat,cbind(scenario=i,design=j,results[[k]]$outputType1))
      }
      k=k+1
    }
  }
  return(list(size=sampleSizeMat,power=powerMat,futility=futilityMat,responseRate=responseRateMat,duration=durationMat,type1=type1Mat))
}

getScenarioSummary=function(isaScenarioList,nISA){
out=array(0,dim=c(nISA,4))
isaScenarioList=isaScenarioList[[3]]
controlRate=isaScenarioList[[2]]$unfavorableEventProb
for(i in seq(3,nISA)){
	out[i,]=(isaScenarioList[[i]]$unfavorableEventProb-isaScenarioList[[2]]$unfavorableEventProb)>=0.0999
}
return(out)
}