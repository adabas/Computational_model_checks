function [data, BIC, iBEST, BEST, pars, NegLL] = subjectAnalysis_v2(subID, datapath, nMod, pbound)

% SUBJECTANALYIS_v2
% This function loads and restructures the subject data. The restructured
% data is fed into a function that determines the fit of the model to the
% data.
%
% INPUTS:
%       subID       : subject id as double
%       datapath    : root directory of the subject folders as char var.
%       nMod        : number of models to test
%       pbound      : lower and upper parameter bounds for each model
%
% OUTPUT:

%% Section 1: Load data

[data, idPartial] = loadSubjectData(subID, datapath);

%% Section 2: Save best fitting parameter and the BIC for each of the models

[BIC, iBEST, BEST, pars, NegLL] = fit_all_v2(data.stimuli, data.rate.binary, idPartial, nMod, pbound, data.stimPresented);

end