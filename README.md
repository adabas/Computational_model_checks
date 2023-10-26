# LSim_simulationChecks

This repository contains the codes for running simulations and checks on the models for the LSim project.

## Script overview
To run the various simulations and checks, run the code in the order of the "Steps". At the starting of each script, you will find commented information about the script.

In brief,
1. *Step1_simulations.m* : provides information about the models of interest. Simulates fake data from the models and plots some model-independent measures.
2. *Step3_parameterRecovery.m* : plots how well the parameters are recovered for the simulated models.
3. *Step4_modelRecovery.m* : sanity checks to validate the models of interest.
4. *Step5_plotData_subjectLevel_v2.m* : to plot subject data. NOTE, this script will NOT run because of the missing participant data file.
5. *Step6_plotData_GroupLevel_v2*: model fit for all participants. Again, will not run due to missing participant data paths.

## Coded on
Matlab version 9.10
