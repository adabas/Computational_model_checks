function [LL, iBEST, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta)

% LIK_ALL_v1
% Function to compute the log likelihood of the data given the parameters
% of a model. Currently, we have three models of interest: the random model
% with one free parameter, the noisy WSLS model with one free paramter and
% the Rescorla Wagner + softmax function model with two free parameters.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pt      : vector containing partial trial numbers
%       par     : 1Xn vector containing parameter values where n refers to
%                 the number of parameters
%
% OUPUT:
%       BIC     : a 1Xn vector containing the BIC values for each of the models.
%       iBEST   : the model number that best accounts for the data
%       BEST    : a 1Xn vector containing information about the degree to
%                 which each of the models best account for the data. Useful
%                 when creating the confusion matrix.
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   May 2020
% =========================================================================

% store the log likelihood values
LL(1) = lik_M1random_v1(a, b);
LL(2) = lik_M2WSLS_v1(a, r, epsilon);
LL(3) = lik_M3RescorlaWagner_v1(a, r, alpha, beta, pt);

% determine the model with the best likelihood value
[M, iBEST] = min(LL);
BEST = LL == M;

end