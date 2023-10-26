function [NegLL, PP, delta, QQ, CK] = lik_M4RWCK(a, r, alpha, beta, alpha_c, beta_c, s)

% LIK_M4RWCK
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices
%       r       : reward
%       alpha   : RW alpha
%       beta    : RW beta
%       alpha_c : CK alpha
%       beta_c  : CK beta
%       s       : stimuli presented trial wise
%
% OUPUT:
%       NegLL   : the negative log likelihood value
%       PP      : trial-wise choice probability matrix
%       delta   : trial-wise prediction error vector
%       QQ      : trial-wise choice value matrix
%       CK      : trial-wise choice kernel matrix
%
% Aroma Dabas
% October 2022
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

% sort presented stimuli as [HR LR]
sSorted = sort(s,2);

% loop over all trial
for t = 1:T
    
    % store value
    QQ(t,:) = q;
    CK(t,:) = k;
    
    % k and q for presented
    q_sub = q(sSorted(t,:));
    k_sub = k(sSorted(t,:));
    
    % compute choice probabilities
    p = M4_softmaxRWCK(q_sub, k_sub, beta, beta_c);
    
    % store choice probability
    PP(t,:) = p;
    
    if isnan(a(t)) || a(t) == 0
        choiceProb(t) = NaN;
        delta(t) = NaN;
    
    else
        % compute choice probability for actual choice
        choiceProb(t) = p(sSorted(t,:) == a(t));

        % update value and choice kernel
        [q(a(t)), k(a(t)), delta(t)] = M4_valueUpdate(alpha, alpha_c, q(a(t)), k(a(t)), r(t));
    end
    
end

% for the last trial
QQ(t+1,:) = q;
CK(t+1,:) = k;

% remove missed trial choice value, kernel and probability
QQ(find(a == 0)+1,:) = [];
CK(find(a == 0)+1,:) = [];
PP(a == 0,:) = NaN;

% compute negative log-likelihood
NegLL = -sum(log(choiceProb(~isnan(choiceProb))));

end