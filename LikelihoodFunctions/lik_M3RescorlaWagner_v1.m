function [NegLL, PP, delta, QQ] = lik_M3RescorlaWagner_v1(a, r, alpha, beta, pt)

% LIK_M3RESCORLAWAGNER_v1
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
%
% Aroma Dabas
% January 2020
% =========================================================================

% store the initial (liking) value of the stimuli
Q = [0.5 0.5];

% number of trials
T = length(a);

% store the evolving probabilities and prediction error
PP = nan(T, 2);
delta = nan(1, T);
QQ = nan(T,2);


% loop over all trial
for t = 1:T
    
    % store the value
    QQ(t,:) = Q; 
    
    % compute choice probabilities
    ev  = exp(beta*Q);
    sev = sum(ev);
    p   = ev / sev;
    
    % store choice probabilities
    PP(t,:) = p;
    
    % compute choice probability for actual choice
    choiceProb(t) = p(a(t));
    
    % update values but not for partial trials
    if ismember(t, pt)
        delta(t) = 0;
    else
        delta(t)   = r(t) - Q(a(t));
    end
    Q(a(t)) = Q(a(t)) + alpha * delta(t);

end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
