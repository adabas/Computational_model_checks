function [BIC, iBEST, BEST] = fit_all_v1(a, r, pt)

% FIT_ALL_v1
% Function to compute the parameter values of each of the models of
% interest that best fit the data. Currently, we have three models of
% interest: the random model with one free parameter, the noisy WSLS model
% with one free paramter and the Rescorla Wagner + softmax function model
% with two free parameters.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pt      : vector containing partial trial numbers
%
% OUPUT:
%       BIC     : a 1X(number of models) vector containing the BIC values
%                 for each of the models.
%       iBEST   : the model number that best accounts for the data
%       BEST    : a 1X(number of models) vector containing information
%                 about the degree to which each of the models best
%                 account for the data. Useful when creating the confusion
%                 matrix.
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   January 2020
% =========================================================================



[~, ~, BIC(1)] = fit_M1random_v1(a);
[~, ~, BIC(2)] = fit_M2WSLS_v1(a, r);
[~, ~, BIC(3)] = fit_M3RescorlaWagner_v1(a, r, pt);
% [~, ~, BIC(4)] = fit_M4CK_v1(a, r);
% [~, ~, BIC(5)] = fit_M5RWCK_v1(a, r);

[M, iBEST] = min(BIC);
BEST = BIC == M;

% why is this step neccesary? Isn't it always 1?
BEST = BEST / sum(BEST);

end