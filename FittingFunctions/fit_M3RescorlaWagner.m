function [Xfit, NegLL, b] = fit_M3RescorlaWagner(a, r, pbound, s)

% FIT_M3RESCORLAWAGNER
% Function to find the parameter values that best fit the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pbound  : parameter bounds [lower; upper];
%       s       : trial wise stimuli presentation
%
% OUPUT:
%       Xfit    : a vector containing the best fitting parameter values
%       NegLL   : the neg loglikelihood value for the best fitting parameter values
%       BIC     : the bayesian information criterion value
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% set fmincon settings
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% create function capturing the RW with softmax function
obFunc = @(x) lik_M3RescorlaWagner(a, r, x(1), x(2), s);

% create vector storing random alpha and beta starting values [alpha beta]
X0 = [rand exprnd(4)];

% store the lower and upper bounds of [alpha beta] parameters
LB = pbound(1,:);
UB = pbound(2,:);

% estimate the parameter fit and the neg log
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% store log likelihood
LL = -NegLL;

% calculate the BIC
b = BIC(length(X0), length(a), NegLL);

end