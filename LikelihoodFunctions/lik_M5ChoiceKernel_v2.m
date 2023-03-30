function [NegLL, PP, k] = lik_M5ChoiceKernel_v2(a, alpha_c, beta_c, s)

% LIK_M3RESCORLAWAGNER_v1
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       alpha_c : parameter alpha value
%       beta_c  : parameter beta value
%       s       : stimuli presented on each trial
%
% OUPUT:
%       NegLL   : the negative log likelihood value
%       PP      : matrix containing choice probabilities at each trial
%       k       : matrix containing choice kernel at each trial
%
% Aroma Dabas
% January 2020
% =========================================================================

% number of trials
T = length(a);

% initialise choice kernel
k = [0 0 0 0];
CK = nan(T, size(k, 2));
PP = nan(T, 2); %size(k, 2));
choiceProb = nan(1,T);

% sort stimuli presentation into [HR LR]
sSorted = sort(s, 2);

% loop over all trial
for t = 1:T
    
    % store the value
    CK(t,:) = k;
    
    % ck for presented
    k_sub = k(sSorted(t,:));

    % compute choice probabilities
    p = M5_softmaxCK(k_sub, beta_c);
    
    % store choice probabilities
    PP(t,:) = p;
    
    if isnan(a(t)) || a(t) == 0
        choiceProb(t) = NaN;
        
    else
        % compute choice probability for actual choice
        choiceProb(t) = p(sSorted(t,:) == a(t));

        % update choice kernel
        k(a(t)) = M5_CKUpdate(alpha_c, k(a(t)));
    end

end

CK(t+1,:) = k;
CK(find(a == 0)+1,:) = [];

PP(a == 0,:) = NaN;

% compute negative log-likelihood
NegLL = -sum(log(choiceProb(~isnan(choiceProb))));

end
