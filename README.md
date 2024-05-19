# Computational analyses: Imagined experiences induce endogenous reinforcement learning

This repository contains the codes for:
* simulating data with the five models of interest. The scripts also run sanity checks on the code.
* computing model fit to the participant data (choices and rewards).

## Script overview
To run each script, check and, if required, update the first section.

In brief,
1. *Step1_simulations.m* : simulates data using each of the five models.
2. *Step2_parameterRecovery.m* : sanity check to test if the code can recover the parameter that was used to generate the data.
3. *Step3_modelRecovery.m* : sanity checks to test if the code can recover the model that was used to generate the data.
4. *Step4_fitModel_subjectLevel.m* : fit models to each of the participant's choices and rewards.
5. *Step5_fitModel_groupLevel*: fit models to all the participants to test model fit at the group level.

## Coded on
Matlab version 9.10
