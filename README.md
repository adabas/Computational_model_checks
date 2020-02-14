# LSim_simulationChecks

This repository contains the codes for running simulations and checks on the models for the LSim project.

## Folder structure
To understand the folder structures in this repository, refer to the "folderStructure.xls" document (you may need to download it).

## Script overview
To run the various simulations and checks, run the code in the order of the "Steps" prepended to the scripts' name. At the starting of each script, you will find commented information about the script.

In brief,
1. *Step1_simulations.m* : simulates fake data from the models and plots some model-independent measures.
2. *Step2_fminconGridSearch.m* : estimates parameters for the RW+softmax function using both fmincon and grid search approaches. It is a sanity check.
3. *Step3a_parameterRecovery.m* : plots how well the parameters are recovered for the RW+softmax function
4. *Step3b_parameterRecovery.m* : plots the bias and variance of the estimated parameters compared to the parameters used for simulating data (simulation parameters). This is done for all the models.
5. *Step4_modelRecovery.m* : sanity checks to validated the models of interest.
6. *Step5_plotData_subjectLevel.m* : to plot subject data. NOTE, this script will NOT run because of the missing participant data file.

## Coded on
Matlab version 9.3
