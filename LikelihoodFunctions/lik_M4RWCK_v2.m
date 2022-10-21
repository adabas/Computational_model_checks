function [NegLL, PP, delta, QQ, CK] = lik_M4RWCK_v2(a, r, alpha, beta, alpha_c, beta_c, pt, stimPres)

% LIK_M4RWCK_v1
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       alpha   : parameter alpha value
%       beta    : parameter beta value
%       pt      : vector containing partial trial numbers
%
% OUPUT:
%       NegLL   : the negative log likelihood value
%       PP      : matrix containing choice probabilities at each trial
%       delta   : vector containing prediction error at each trial
%       QQ      : matrix containing choice values at each trial
%       CK      : matrix containing choice kernel value update
%
% Aroma Dabas
% April 2020
% =========================================================================

% number of trials
T = length(a);

% initialise values
q  = [0.5 0.5 0.5 0.5]; % initial expected reward values
k = [0 0 0 0];     % initialise choice kernel
PP = nan(T, 2);%PP = nan(T, size(q, 2)); % probabilities
QQ = nan(T, size(q, 2)); % value update
CK = nan(T, size(q, 2)); % choice kernal update
delta = nan(T, 1); % prediction error

% loop over all trial
for t = 1:T
    
    % store value
    QQ(t,:) = q;
    CK(t,:) = k;
    
    % k and q for presented
    q_sub = q(stimPres(t,:));
    k_sub = k(stimPres(t,:));
    
    % compute choice probabilities
    p = M4_softmaxRWCK(q_sub, k_sub, beta, beta_c);
    
    % store choice probability
    PP(t,:) = p;
                
    % compute choice probability for actual choice
    choiceProb(t) = p(stimPres(t,:) == a(t));%p(a(t));
    
    % update value and choice kernel
    [q(a(t)), k(a(t)), delta(t)] = M4_valueUpdate(alpha, alpha_c, q(a(t)), k(a(t)), r(t), t, pt);
    
end

% for the last trial
QQ(t,:) = q;
CK(t,:) = k;

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));

end