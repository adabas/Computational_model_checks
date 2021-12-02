function [bic, iBEST, BEST, pars, NegLL] = fit_all_v1(a, r, pt, nMod, pbound, keycode)

% FIT_ALL_v1 Function to compute the parameter values of each of the models
% of interest that best fit the data. Currently, we have five models of
% interest: the random model with one free parameter, the noisy WSLS model
% with one free paramter, the Rescorla Wagner + softmax function model with
% two free parameters, RW with choice kernel model with four parameters and
% a simple choice kernel model with two free parameters.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pt      : vector containing partial trial numbers
%       nMod    : number of models
%       pbound  : parameter bounds [lower; upper];
%       keycode : keypresses
%
% OUPUT:
%       bic     : a 1X(number of models) vector containing the BIC values
%                 for each of the models.
%       iBEST   : the model number that best accounts for the data
%       BEST    : a 1X(number of models) vector containing information
%                 about the degree to which each of the models best
%                 account for the data. Useful when creating the confusion
%                 matrix.
%       pars    : a 5X5 matrix containing the best fitting parameters for
%                 each model.
%       LL      : the negative log likelihood value of the best fitting
%                 parameter.
%
% Aroma Dabas [dabas@cbs.mpg.de]
% January 2020
% =========================================================================

% remove all NaN inputs
if any(isnan(a))
    [id,~] = find(isnan(a));
    a(id,:) = [];
    r(id,:) = [];
end

% iterate to find the parameter values that best fit the data
for iter = 1:10  % run this for 20 iterations
       [x1(iter), l1(iter), b1(iter)] = fit_M1random_v1(keycode, pbound(:,1));
        [x2(iter), l2(iter), b2(iter)] = fit_M2WSLS_v1(a, r, pbound(:,2));
        [x3(iter,:), l3(iter), b3(iter)] = fit_M3RescorlaWagner_v1(a, r, pt, pbound(:,3:4));
        [x4(iter,:), l4(iter), b4(iter)] = fit_M4RWCK_v1(a, r, pt, pbound(:, 3:6));
        [x5(iter,:), l5(iter), b5(iter)] = fit_M5ChoiceKernel_v1(a, r, pt, pbound(:, 5:6));
end

% determine the best fitting parameter for each model
[NegLL(1),i(1)]=min(l1(:));
[NegLL(2),i(2)]=min(l2(:));
[NegLL(3),i(3)]=min(l3(:));
[NegLL(4),i(4)]=min(l4(:));
[NegLL(5),i(5)]=min(l5(:));

% select the BIC for best fitting parameter
bic(1) = b1(i(1));
bic(2) = b2(i(2));
bic(3) = b3(i(3));
bic(4) = b4(i(4));
bic(5) = b5(i(5));

% store the best fitting parameter
pars = nan(nMod);
pars(1,1) = x1(i(1));
pars(2,1) = x2(i(2));
pars(3,1:2) = x3(i(3),:);
pars(4,1:4) = x4(i(4),:);
pars(5,1:2) = x5(i(5),:);

% find the model that best fits the data
[M, iBEST] = min(bic); %min(NegLL);
BEST = bic == M; %NegLL == M;
BEST = BEST / sum(BEST);

end