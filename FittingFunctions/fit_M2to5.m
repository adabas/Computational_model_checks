function [bic, iBEST, BEST, pars, BESTnll] = fit_M2to5(a, r, nMod, pbound, s)

% FIT_M2to5
% Function to compute the parameter values of all models of interest,
% except random responding model, that best fit the data. Currently, we
% have five models of interest: the random model with one free parameter,
% the noisy WSLS model with one free paramter, the Rescorla Wagner +
% softmax function model with two free parameters, RW with choice kernel
% model with four parameters and a simple choice kernel model with two free
% parameters.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       nMod    : number of models
%       pbound  : parameter bounds [lower; upper];
%       s       : stimuli presented at each trial
%
% OUTPUT:
%       bic     : vector of BIC values for each model
%       iBEST   : number of the model that best fits the data
%       BEST    : vector with winning model's index
%       pars    : matrix of the best fitting parameters/model
%       LL      : the negative log-likelihood value of the best fitting
%                 model
%       BESTnll : vector with Neg. LL winning model's index
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% iterate to find the parameter values that best fit the data
for iter = 1:10
    [x2(iter), l2(iter), b2(iter)] = fit_M2WSLS(a, r, pbound(:,2), s);
    [x3(iter,:), l3(iter), b3(iter)] = fit_M3RescorlaWagner(a, r, pbound(:,3:4), s);
    [x4(iter,:), l4(iter), b4(iter)] = fit_M4RWCK(a, r, pbound(:, 3:6), s);
    [x5(iter,:), l5(iter), b5(iter)] = fit_M5ChoiceKernel(a, r, pbound(:, 5:6), s);
end

% determine the best fitting parameter for each model
[NegLL(1),~]=min(l2(:));
[NegLL(2),~]=min(l3(:));
[NegLL(3),~]=min(l4(:));
[NegLL(4),~]=min(l5(:));

% select the BIC for best fitting parameter
[bic(1),i(1)]=min(b2(:));
[bic(2),i(2)]=min(b3(:));
[bic(3),i(3)]=min(b4(:));
[bic(4),i(4)]=min(b5(:));

% store the best BIC fitting parameter
pars = nan(nMod);
pars(1,1) = x2(i(1));
pars(2,1:2) = x3(i(2),:);
pars(3,1:4) = x4(i(3),:);
pars(4,1:2) = x5(i(4),:);

% find the model that best fits the data using the BIC
[M, iBEST] = min(bic);
BEST = bic == M;
BEST = BEST / sum(BEST);

% find the model that best fits the data using the NegLL
[Mnll, iBESTnll] = min(NegLL);
BESTnll = NegLL == Mnll;
BESTnll = BESTnll / sum(BESTnll);

end