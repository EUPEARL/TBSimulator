######
### Programmed by: Tobias  Mielke
### These functions are specific to TB for the definition of response scenarios on endpoints
### Version: Draft (run)
######


isaDesignParameters=list()

###
# Parameter 1: Probability of culture conversion
# Parameter 2: Probability of conversion at follow up (if responder)
# Parameter 3: Probability of non-favorable event
# Parameter 4: Probability of non-favorable event (if non-responder)
isaSimulationScenarios=list()
worstCase=c(0.999,0.8,0.05)
bestCase=c(0.999,0.9,0.05)
shapeParameter=c(1.65)
followUp=c(12)

### Bad Good
relativeEffect=0
controlScenario=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)


worstCase=c(0.999,0.8,0.15)
bestCase=c(0.999,0.9,0.05)

### Bad Base
relativeEffect=rep(0,4)
baseEffect0=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Ok Base
relativeEffect=rep(0.5,4)
baseEffect05=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Good Base
relativeEffect=rep(1,4)
baseEffect1=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Bad Linear
relativeEffect=rep(0,4)
linearEffect0=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Ok Linear
relativeEffect=seq(0,0.5,length=4)
linearEffect05=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Good Linear
relativeEffect=seq(0,1,length=4)
linearEffect1=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Bad Expo
relativeEffect=rep(0,4)
expoEffect0=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Ok Expo
relativeEffect=c(0,0,0,0.5)
expoEffect05=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

### Good Expo
relativeEffect=c(0,0,0,1)
expoEffect1=getEffectScenario(relativeEffect,worstCase,bestCase,shapeParameter,followUp)

