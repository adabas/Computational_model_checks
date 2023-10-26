function [bic, iBEST, BEST, pars, NegLL] = fit_all_models(a, r, nMod, pbound, s)

% FIT_ALL_MODELS
% Function to compute the parameter values of each of the models
% of interest that best fit the data. Currently, we have five models of
% interest: the random model with one free parameter, the noisy WSLS model
% with one free paramter, the Rescorla Wagner + softmax function model with
% two free parameters, RW with choice kernel model with four parameters and
% a simple choice kernel model with two free parameters.
%
% Difference from v1: Updated to fit data to the four stimuli, instead of
% the two categories.
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
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% iterate to find the parameter values that best fit the data
for iter = 1:10
   [x1(iter), l1(iter), b1(iter)] = fit_M1random(a, pbound(:,1), s);
    [x2(iter), l2(iter), b2(iter)] = fit_M2WSLS(a, r, pbound(:,2), s);
    [x3(iter,:), l3(iter), b3(iter)] = fit_M3RescorlaWagner(a, r, pbound(:,3:4), s);
    [x4(iter,:), l4(iter), b4(iter)] = fit_M4RWCK(a, r, pbound(:, 3:6), s);
    [x5(iter,:), l5(iter), b5(iter)] = fit_M5ChoiceKernel(a, r, pbound(:, 5:6), s);
end

% determine the best fitting parameter for each model
[NegLL(1),~]=min(l1(:));
[NegLL(2),~]=min(l2(:));
[NegLL(3),~]=min(l3(:));
[NegLL(4),~]=min(l4(:));
[NegLL(5),~]=min(l5(:));

% select the BIC for best fitting parameter
[bic(1),i(1)]=min(b1(:));
[bic(2),i(2)]=min(b2(:));
[bic(3),i(3)]=min(b3(:));
[bic(4),i(4)]=min(b4(:));
[bic(5),i(5)]=min(b5(:));

% store the best BIC fitting parameter
pars = nan(nMod);
pars(1,1) = x1(i(1));
pars(2,1) = x2(i(2));
pars(3,1:2) = x3(i(3),:);
pars(4,1:4) = x4(i(4),:);
pars(5,1:2) = x5(i(5),:);

% find the model that best fits the data
[M, iBEST] = min(bic);
BEST = bic == M;
BEST = BEST / sum(BEST);

end