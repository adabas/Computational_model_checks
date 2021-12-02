function [NegLL, PP, k] = lik_M5ChoiceKernel_v1(a, alpha_c, beta_c)

% LIK_M3RESCORLAWAGNER_v1
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       alpha_c : parameter alpha value
%       beta_c  : parameter beta value
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
k = [0 0];
CK = nan(T, size(k, 2));
PP = nan(T, size(k, 2));
choiceProb = nan(1,T);

% loop over all trial
for t = 1:T
    
    % store the value
    CK(t,:) = k;

    % compute choice probabilities
    p = M5_softmaxCK(k, beta_c);
    
    % store choice probabilities
    PP(t,:) = p;
    
    % compute choice probability for actual choice
    choiceProb(t) = p(a(t));
    
    % update choice kernel
    %k(a(t)) = M5_CKUpdate(alpha_c, k(a(t)));
    k = (1-alpha_c) * k;
    k(a(t)) = k(a(t)) + alpha_c * 1;

end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
