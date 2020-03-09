function [Xfit, LL, BIC] = fit_M3RescorlaWagner_v1(a, r, pt)

% FIT_M3RESCORLAWAGNER_v1
% Function to find the parameter values that best fit the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pt      : vector containing partial trial numbers
%
% OUPUT:
%       Xfit    : a vector containing the best fitting parameter values
%       LL      : the negative loglikelihood value for the best fitting
%                 parameter values
%       BIC     : the bayesian information criterion value
%
% Aroma Dabas [dabas@cbs.mpg.de]
% January 2020
% =========================================================================

% set fmincon settings
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% create function capturing the RW with softmax function
obFunc = @(x) lik_M3RescorlaWagner_v1(a, r, x(1), x(2), pt);

% create vector storing random alpha and beta starting values [alpha beta]
X0 = [rand exprnd(1)];

% store the lower and upper bounds of [alpha beta] parameters
LB = [0 0];
UB = [1 20];

% estimate the parameter fit and the neg log
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% store the log likelihood as the negative log likelihood
LL = -NegLL;

% calculate the BIC
BIC = length(X0) * log(length(a)) + 2*NegLL;

end