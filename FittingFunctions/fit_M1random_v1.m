function [Xfit, LL, b] = fit_M1random_v1(a, pbound)

% FIT_M1RANDOM_v1
% Function to compute the parameter values that best fit the data.
%
% INPUT:
%       a       : choices vector
%       pbound  : parameter bound [lower; upper];
%
% OUPUT:
%       Xfit    : a vector containing the best fitting parameter value
%       LL      : the loglikelihood value for the best fitting parameter value
%       b       : the bayesian information criterion value
%
% Aroma Dabas [dabas@cbs.mpg.de]
% January 2020
% =========================================================================

% set fmincon settings
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% Find the best fitting parameters and the negative loglikelihood, while
% fitting a function (obFunc) to the data (a is the choice and r is the
% reward).

% specify the function
obFunc = @(x) lik_M1random_v1(a, x);

% random starting point
X0 = rand;

% bounds
LB = pbound(1);
UB = pbound(2);

% run fmincon
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% save the values as negative log
LL = -NegLL;

% calculate the BIC
b = BIC(length(X0), length(a), NegLL);

end