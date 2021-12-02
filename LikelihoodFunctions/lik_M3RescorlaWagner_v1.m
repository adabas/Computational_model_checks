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

% number of trials
T = length(a);

% initialise values
q = [0.5 0.5];      % initial value of both stimuli
PP = nan(T,2);      % probability of selecting either of the stimuli
delta = nan(1,T);   % prediction error
QQ = nan(T, size(q, 2));    % trial wise updated value
choiceProb = nan(1,T);      % trial wise choice probability

% loop over all trial
for t = 1:T
    
    % store the value
    QQ(t,:) = q; 
    
    % compute choice probabilities
    p = M3_softmaxFunction(q, beta);
    
    % store choice probabilities
    PP(t,:) = p;
    
    if isnan(a(t))
        if t == 1   % missed first trial
            q(t) = 0.5;
            delta(t) = 0; %NaN; better 0 - no diff. b/w exp. and actual rew
        else
            choiceProb(t) = choiceProb(t-1);
            delta(t) = 0;
        end
        
    else
        % compute choice probability for actual choice
        choiceProb(t) = p(a(t));

        % value update
        [q(a(t)), delta(t)] = M3_valueUpdate(alpha, q(a(t)), r(t), t, pt);
    end
end

% for last trial
QQ(t,:) = q; 

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
