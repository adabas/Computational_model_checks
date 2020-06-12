function [BIC, iBEST, BEST, pars, LL] = fit_all_v1(a, r, pt, nMod, pbound)

% FIT_ALL_v1
% Function to compute the parameter values of each of the models of
% interest that best fit the data. Currently, we have three models of
% interest: the random model with one free parameter, the noisy WSLS model
% with one free paramter and the Rescorla Wagner + softmax function model
% with two free parameters.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       pt      : vector containing partial trial numbers
%       nMod    : number of models
%       pbound  : parameter bounds [lower; upper];
%
% OUPUT:
%       BIC     : a 1X(number of models) vector containing the BIC values
%                 for each of the models.
%       iBEST   : the model number that best accounts for the data
%       BEST    : a 1X(number of models) vector containing information
%                 about the degree to which each of the models best
%                 account for the data. Useful when creating the confusion
%                 matrix.
%       pars    : a 4X4 matrix containing the best fitting parameters for
%                 each model.
%       LL      : the negative log likelihood value of the best fitting
%                 parameter.
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   January 2020
% =========================================================================

% iterate to find the best fitting parameter
for iter = 1:10
    for m = 1:nMod
        if m == 1
            [x1(iter), l1(iter), b1(iter)] = fit_M1random_v1(a, pbound(:,m));
        elseif m == 2
            [x2(iter), l2(iter), b2(iter)] = fit_M2WSLS_v1(a, r, pbound(:,m));
        elseif m == 3
            [x3(iter,:), l3(iter), b3(iter)] = fit_M3RescorlaWagner_v1(a, r, pt, pbound(:,m:m+1));
        elseif m == 4
%           [x4(iter,:), l4(iter), b4(iter)] = fit_M4RWCK_v1(a, r, pt, pbounds(:, m+1:m+4));
        end
    end
end

% determine the best fitting parameter for each model
[LL(1),i(1)]=min(l1(:));
[LL(2),i(2)]=min(l2(:));
[LL(3),i(3)]=min(l3(:));
% [LL(4),i(4)]=min(l4(:));

% select the BIC for best fitting parameter
BIC(1) = b1(i(1));
BIC(2) = b2(i(2));
BIC(3) = b3(i(3));
% BIC(4) = b4(i(4));

% store the best fitting parameter
pars = nan(nMod);
pars(1,1) = x1(i(1));
pars(2,1) = x2(i(2));
pars(3,1:2) = x3(i(3),:);
% pars(4,1:4) = x4(i(4),:);

% Note: clean script from line 49 till 65

% find the model that best fits the data
[M, iBEST] = min(BIC);
BEST = BIC == M;
BEST = BEST / sum(BEST);

end