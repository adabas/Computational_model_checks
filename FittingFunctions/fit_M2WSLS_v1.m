function [Xfit, LL, b] = fit_M2WSLS_v1(a, r)

% FIT_M2WSLS_v1
% Function to compute the parameter values that best fit the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%
% OUPUT:
%       Xfit    : a vector containing the best fitting parameter value
%       LL      : the loglikelihood value for the best fitting parameter value
%       BIC     : the bayesian information criterion value
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   January 2020
% =========================================================================

% set fmincon settings
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% Find the best fitting parameters and the negative loglikelihood, while
% fitting a function (obFunc) to the data (a is the choice and r is the
% reward).

% specify the function
obFunc = @(x) lik_M2WSLS_v1(a, r, x);

% random starting point
X0 = rand;

% state the upper and lower bounds
LB = 0;
UB = 1;

% estimate the parameters
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% store the negative log
LL = -NegLL;

% calculate the BIC
b = BIC(length(X0), length(a), NegLL);
%BIC = length(X0) * log(length(a)) + 2*NegLL;

end