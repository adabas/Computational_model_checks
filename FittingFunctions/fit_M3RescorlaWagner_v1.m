function [Xfit, LL, BIC] = fit_M3RescorlaWagner_v1(a, r, pt)

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