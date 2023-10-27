function [data, BIC, iBEST, BEST, pars, NegLL] = subjectAnalysis(subID, datapath, nMod, pbound)

% ----------------------------------------------------------------------------------------
% SUBJECTANALYIS
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
%		data		: structure of participant's trail-wise choices and feedback
%		BIC			: vector of BIC values from each model fit
%		iBEST 		: the number corresponding to the model best fitting the data
%		BEST		: vector of 0s and 1 (index of the model best fitting the data)
%		pars		: matrix containing best fitting parameters for each model
%		NegLL		: vector of negative log likelihoods from model fits
%
% Aroma Dabas [dabas@cbs.mpg.de]
% ----------------------------------------------------------------------------------------

%% Section 1: Load data

[data, idPartial] = loadSubjectData(subID, datapath);

%% Section 2: Save best fitting parameter and the BIC for each of the models

[BIC, iBEST, BEST, pars, NegLL] = fit_all_models(data.stimuli, data.rate.binary, nMod, pbound, data.stimPresented);

end