##### COPYRIGHT #############################################################################################################
#
# Copyright (C) 2023 JANSSEN RESEARCH & DEVELOPMENT, LLC
# This package is governed by GNU General Public License V3 with additional terms. The precise license terms are located in the files
# LICENSE and GPL.
#
#############################################################################################################################.

#############################################################################################################################.
#   Description                                                                                                         
#   This file contains functions which build the API to the customized code to conduct analyses, get decisions and implement those for the ISAs.
#   Developer(s): T. Mielke PhD                                                                               			
#############################################################################################################################.

# conductISAAnalysis
# Input: 
#	time: Time in the platform
#	interim: Integer indicating which interim analysis is being conducted
#	isaDesign: the ISA design list containing parameters of relevance for the specific ISA
#	data:		analysis data for the specific ISA only (containing already control data)
# Output:
#	list with the first element being the interim analysis and the second element being the analysisResults as passed by the customAnalysisFunction
conductISAAnalysis=function(time,interim,isaDesign,data){
	nEndpoints=isaDesign$nEndpoints
	analysisResults=list()
	analysisResults[[1]]=customAnalysisFunctions(time,data,interim,isaDesign)
	
	return(list(interim=interim,
					analysisResults=analysisResults))
}

# analysisWrapper
# Input: 
#	time: Time in the platform
#	isaDesignList: The complete list of all ISA designs which have already been initiated already.
#	isaData: a list with analysis data for each ISA (containing already control data)
#	conductAnalysis: a vector indicating which analysis is to be conducted per ISA. If "0", no analysis will be conducted.
# Output:
#	a list with analysisResults for each ISA. Each of the analysis results consists of the interim and the actual results from the analysis.
analysisWrapper=function(time,isaDesignList,isaData,conductAnalysis){
	analysisResults=list()
	if(any(conductAnalysis>0)){
		for(i in seq(1,length(conductAnalysis))){
			if(conductAnalysis[i]>0){
				analysisResults[[i]]=conductISAAnalysis(time,conductAnalysis[i],isaDesignList[[2+i]],isaData[[i]])
			}
		}
	} else {
		return(NULL)
	}
	return(analysisResults)
}

# getDecisions
# Input: 
#	conductAnalysis: a vector indicating which analysis is to be conducted per ISA. If "0", no analysis will be conducted.
#	isaDesignList: The complete list of all ISA designs which have already been initiated already.
#	analysisResults: a list with analysisResults for each ISA. Each of the analysis results consists of the interim and the actual results from the analysis.
# Output:
#	decisionResults: a list with decisions for each ISA with an analysis, as specified in getDecisionISA
getDecisions=function(conductAnalysis,isaDesignList,analysisResults){
	decisionResults=list()
	if(any(conductAnalysis>0)){
		for(i in seq(1,length(conductAnalysis))){
			if(conductAnalysis[i]>0){
				#print(c("ISA:",i))
				decisionResults[[i]]=getDecisionISA(isaDesignList[[i+2]],analysisResults[[i]])
			}
		}
		return(decisionResults)
	} else {
		return(NULL)
	}
}

# getDecisionsISA
# Input: 
#	isaDesign: the ISA design list containing parameters of relevance for the specific ISA
#	analysisResults: a list with the counter for the interim analysis and the analysisResults for the respective ISA.
# Output:
#	decisions: a list with decisions for the ISA, containing at minimum: summary, sampleSize, cutoffs, population (not used), arms (not used), allocProb
getDecisionISA=function(isaDesign,analysisResults){
	#print(analysisResults)
	interim=analysisResults[[1]]
	analysisResults=analysisResults[[2]]
	decisions=getDecisionCustom(interim,isaDesign,analysisResults)
}

# implementDecisions
# Input: 
#	time: time in platform
#	conductAnalysis: a vector indicating which analysis is to be conducted per ISA. If "0", no analysis will be conducted.
#	isaDesignList: The complete list of all ISA designs which have already been initiated already.
#	analysisResults: a list with analysisResults for each ISA. Each of the analysis results consists of the interim and the actual results from the analysis.
#	decisions: a list with decisions for each ISA
# Output:
#	isaDesignList: the updated isaDesignList, which contains new design rules and data based on the interim analysis.
implementDecisions=function(time,conductAnalysis,isaDesignList,analysisResults,decisions){
	if(any(conductAnalysis>0)){
			for(i in seq(1,length(conductAnalysis))){
				if(conductAnalysis[i]>0){
				isaDesignList[[i+2]]=implementDecisionsISA(time,isaDesignList[[i+2]],analysisResults[[i]],decisions[[i]])
			}
		}
	}
	return(isaDesignList)
}

# implementDecisionsISA
# Input: 
#	time: time in platform
#	isaDesign: the specific ISA design
#	analysisResults: the analysisResults for the specific ISA (first element interim counter, following element analysis result)
#	decision: a list with decisions for the specific ISA
# Output:
#	isaDesign: the updated isaDesign, which contains new design rules and data based on the interim analysis.
implementDecisionsISA=function(time,isaDesign,analysisResults,decision){
	isaDesign$stage=isaDesign$stage+1
	isaDesign$interimPassed[1:analysisResults[[1]]]=1
	isaDesign$analysisTime[isaDesign$stage]=time
	isaDesign$targetSize=decision$sampleSize
	isaDesign$interimCutoff=decision$cutoffs
		
	isaDesign$allocProb=decision$allocProb
	### NEED TO ADJUST RANDOMIZATION LIST?

	if(decision$summary[2]==0 | !any(isaDesign$currentSize<isaDesign$targetSize)){
		isaDesign$activeISA=FALSE
	} else {
		isaDesign$activeISA=TRUE
	}
	if((analysisResults[[1]]==length(isaDesign$interimPassed)&isaDesign$interimTrigger!="perpetual"&isaDesign$interimTrigger!="customTrigger")|
		decision$summary[length(decision$summary)]==1){
		isaDesign$activeISA=FALSE
		isaDesign$completed=TRUE
	}
	
	if(isaDesign$interimTrigger=="perpetual"|isaDesign$interimTrigger=="customTrigger"){
		if(isaDesign$stage==2){
			isaDesign$decisions[analysisResults[[1]],]=decision$summary
		} else {
			isaDesign$decisions=rbind(isaDesign$decisions,decision$summary)
		}
	} else {
		isaDesign$decisions[analysisResults[[1]],]=decision$summary
	}
	return(isaDesign)
}