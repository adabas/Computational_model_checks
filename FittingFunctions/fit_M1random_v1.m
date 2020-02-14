function [Xfit, LL, BIC] = fit_M1random_v1(a)

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
LB = 0;
UB = 1;

% run fmincon
[Xfit, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

% save the values as negative log
LL = -NegLL;

% calculate the BIC
BIC = length(X0) * log(length(a)) + 2*NegLL;

end