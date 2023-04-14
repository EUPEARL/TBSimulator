# TBSimulator
A generic simulator has been developed to support design discussions between EU-PEARL WP2 (Methodology) and WP5 (Tuberculosis). The concept of working with generic code was considered to enable larger flexibility in the iterative platform design discussion process. The resulting code consists of two blocks:
- Generic project agnostic code conducting standard tasks, which are required in many platform trials, such as patient recruitment, data generation, monitoring for interim analyses and implementing decisions resulting from interim analyses.
- Project specific code for conducting tasks which are specific to the project, such as parametric time-to-event modelling for Tuberculosis, definition of effect assumptions and specific ISA-designs.

The generic code has been stored in R-scripts, which are separated by functionality:
- allocationFunctions.R: supporting allocation of study participants to the ISAs and within ISAs
- dataGenerationFunctions.R: generation of patient level data
- dataProcessingFunctions.R: processing of data for analyses
- monitoringFunctions.R: monitoring for interim analyses
- simulationFunctions.R: running a single platform trial simulation, including updates of ISAs
- mainSimulationWrapper.R: wrapper for conducting multiple simulations
- genericAnalysisFunctions.r: functions for analyses, decision making and decision implementation - API towards "customFunctions.r"
- supportTools.r: small supporting functions simple operations.
- sourceLibrary: code sourcing all other scripts.

The project specific code covers:
- designSetFunctions: Definition of the ISA and platform designs and simulation settings
- isaScenarios: Defintition of response scenarios for across ISAs
- responseScenarios: Definition of the actual available set of response scenarios
- customFunctions: The API to genericAnalysisFunctions.r, containing code to analyze the data and to make decision
- simulationResultCompiling: Code to post-process outcomes from the platform simulation for calculation of operating characteristics
- sourrceProject: code sourcing all project specific scripts
- simulationSetup: code which executes the simulation

The simulator has been used for simulations on allocation rules, as well as for the simulations of the WP5 Tuberculosis platform design with the aim of regimen shortening. The code is very flexible, but also complex to understand and adjust. The code generates patient level data in the long data format, which may result in relatively large data sets over a platform trial, if many visits are being conducted per patient. This may result in relatively slow computational performance.

Note that the code is using the packages randomizeR (GPL-3), data.table (MPL-2.0) and mvtnorm (GPL-2). The use of data.table has been idenfitied as a potential issue for further optimization, as data.table may result in excessive memory usage.

The code has not been formally validated. Results have been informally assessed for simplistic scenarios against standard clinical trial simulations. Results have been informally cross-validated against SIMPLE for a test case developed on the basis of the EU-PEARL MDD-platform. 
