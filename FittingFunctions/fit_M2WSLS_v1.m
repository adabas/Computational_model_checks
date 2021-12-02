function [Xfit, NegLL, b] = fit_M2WSLS_v1(a, r, pbound)

% FIT_M2WSLS_v1
% Function to compute the parameter values that best fit the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pbound  : parameter bounds [lower; upper];
%
% OUPUT:
%       Xfit    : a vector containing the best fitting parameter value
%       NegLL   : the neg. loglikelihood value for the best fitting parameter value
%       BIC     : the bayesian information criterion value
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   January 2020
% =========================================================================

% set fmincon settings
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% specify the function
obFunc = @(x) lik_M2WSLS_v1(a, r, x);

% random starting point
X0 = rand;

% state the upper and lower bounds
LB = pbound(1);
UB = pbound(2);

% estimate the parameters
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% store the log
LL = -NegLL;

% calculate the BIC
b = BIC(length(X0), length(a), NegLL);

end