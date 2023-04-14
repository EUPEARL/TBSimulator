######
### Programmed by: Tobias  Mielke
### These functions trigger a single simulation or a round of multiple simulations.
### Version: Draft (run)
######
				

rm(list=ls())

output="output"
pathLib="~/genericSimulator/"
source(paste(pathLib,"sourceLibrary.r",sep=""))

pathProject="~/projectCode/"
source(paste(pathProject,"sourceProject.r",sep=""))

#isaDesignParameters=list()
  
### Allocation rules
  r0SQ=c("perActiveArm","SR",0.15,1)
     
  ###
  # allocProb = Prior probability of being non-inferior to control
  # how to improve?
  # Use predictive probability of study success to define allocation ratio instead of post probability?
  # ... predictive probability: more data available => Predictive probability increases
  									
  longCohort=getDesign(nArms=1,sizePerArm=80,trtDuration=26,
  	priorAllocProb=1,postProbCutoff=c(0,0,-1),ia=2,activeISA=FALSE,name="long_cohort",scenario=1,dropOutRate=0.01,allocationRule="constant")
  
  controlISA=getDesign(nArms=1,sizePerArm=80,trtDuration=26,
  	priorAllocProb=1,postProbCutoff=c(0,0,-1),ia=2,activeISA=TRUE,name="control",scenario=1,dropOutRate=0.01,allocationRule="constant")
  
  activeISA_Base=getDesign(nArms=4,sizePerArm=80,trtDuration=c(12,16,20,26),
  	priorAllocProb=rep(0.25,4),postProbCutoff=c(0,0,0,-1),ia=2,activeISA=TRUE,name="active",scenario=1,dropOutRate=0.01,allocationRule="constant")	
	
  activeISA_2IA=getDesign(nArms=4,sizePerArm=80,trtDuration=c(12,16,20,26),
  	priorAllocProb=rep(0.25,4),postProbCutoff=c(0,0,0.1,0.4),ia=2,activeISA=TRUE,name="active",scenario=1,dropOutRate=0.01,allocationRule="constant")	
   	
	isaList=list(longCohort,controlISA,activeISA_Base,activeISA_2IA)
	isaList_Base=defineISAList(isaList,c(1,2,rep(3,8)))
	isaList_2IA=defineISAList(isaList,c(1,2,rep(4,8)))
	
	### Computational settings
	recruitmentSettings=list(list(
	recruitmentParameters=c(0,0,40/4.25,0),recruitmentType="separate"))
		
	computationalSettings=list(randomizationParameter=2,
							blockSize=20,
							useRandomizationList=FALSE,
							simulationSettings=c(104,1,1000,500),
							storeAllocRates=FALSE)

	baselineSimulationParameters=list(
		nCovariates=1,
		covariateType=c("Uniform"),
		covariateNames=c("Age"),
		pars3=c(20,60))
	
	simParams=getSimulationParameters(recruitmentSettings=recruitmentSettings,
				baselineSimulationParameters=baselineSimulationParameters,
				computationalSettings=computationalSettings,
				endpointSettings=list(endpoints=1,endpoint=c("customized")))

### Design List
  platformDesign_Base=getPlatformDesign(allocProb=c(1e-10,1,3)/sum(c(1e-10,1,3)),
  				allocRule=r0SQ,
  				entryTime=c(0,0,seq(0,7*8,by=8)),
				isaOrder=TRUE,isaProb=NULL,randomizer="Target",isaList=isaList_Base)
				  
  platformDesign_r0SQ_2IA=getPlatformDesign(allocProb=c(1e-10,1,3)/sum(c(1e-10,1,3)),
  				allocRule=r0SQ,
				entryTime=c(0,0,seq(0,7*8,by=8)),
  				isaOrder=TRUE,isaProb=NULL,randomizer="Target",isaList=isaList_2IA)  				
  				
  designParams=list(platformDesign_Base,
  	platformDesign_r0SQ_2IA)
  
  
  isaListTMP=list(isaList_Base,isaList_2IA)
  				
  outPath=paste(pathProject,output,"/",sep="")
  
  scenarioName=c("allActive","random")
  designName=c("Base","2IAr0SQ")
  
  ### Specify path to store output from simulations:
  dir.create(outPath)  
	for(j in seq(1,2)){
			dir.create(paste(outPath,designName[j],sep=""))
		for(i in seq(1,2)){  
			dir.create(paste(outPath,designName[j],"/",scenarioName[i],sep=""))
		}
	}
	
  seed=1
  simScen=c(1,2)
  simDes=c(1,2)
  nSim=5	  
  ### For single platform trial simulation
	result=simulateSingleTrial(
				nSim=1,seed=seed,
				platformSimulationParameters=simParams[[1]],
				isaSimulationScenarios=isaScenarioList[[simScen[1]]],
				platformDesign=designParams[[simDes[1]]],
				entryTime=c(0,0,seq(0,7*8,by=8)),
				isaList=isaListTMP[[simDes[1]]])
	
	### For multiple simulations
	results=list()
			kk=1
			for(i in simScen){  
				for(j in simDes){
					print(c("Scenario:",i,"Design:",j))
					 test2IA=simulatePlatform(nSim=nSim,
						  seed=seed,
						  platformSimulationParameters=simParams[[1]],
						  isaSimulationScenarios=isaScenarioList[[i]],
						  platformDesign=designParams[[j]],
							#entryTime=c(0,0,0),
							entryTime=c(0,0,seq(0,7*8,by=8)),
							isaList=isaListTMP[[j]],
							outPath=paste(outPath,designName[j],"/",scenarioName[i],"/",sep=""),
							savingGrid=100)

					results[[kk]]=getDesignCharacteristics(results=NULL,
						path=paste(outPath,designName[j],"/",scenarioName[i],"/",sep=""),
						scenario=paste("summaryOut",seed,sep=""),
						isaList=isaListTMP[[j]])
					kk=kk+1
					}
			}
