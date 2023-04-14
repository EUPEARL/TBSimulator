######
### Programmed by: Tobias  Mielke
### These functions are specific to TB for the definition of the simulation scenarios
### Version: Draft (run)
######

scenarioPool=list(baseEffect0,baseEffect05,baseEffect1,linearEffect0,linearEffect05,linearEffect1,expoEffect0,expoEffect05,expoEffect1)

scenario1=list(randomScenario=FALSE,scenarioProb=NA,scenarioList=list(controlScenario,controlScenario,baseEffect1,baseEffect1,baseEffect1,baseEffect1,baseEffect1,baseEffect1,baseEffect1,baseEffect1))
scenario2=list(randomScenario=TRUE,scenarioProb=c(0,0,rep(1/9,9)),scenarioList=list(controlScenario,controlScenario,baseEffect0,baseEffect05,baseEffect1,linearEffect0,linearEffect05,linearEffect1,expoEffect0,expoEffect05,expoEffect1))

isaScenarioList=list(
scenario1=scenario1,
scenario2=scenario2)
