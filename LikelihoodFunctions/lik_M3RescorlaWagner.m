function [NegLL, PP, delta, QQ, choiceProb] = lik_M3RescorlaWagner(a, r, alpha, beta, s)

% LIK_M3RESCORLAWAGNER
% Function to compute the negative log-likelihood values for fitting the model to the data.
% This script is for stimuli instead of category options.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       alpha   : parameter alpha value
%       beta    : parameter beta value
%       s       : stimuli presented at the trial
%
% OUPUT:
%       NegLL   : the negative log likelihood value
%       PP      : matrix containing choice probabilities at each trial
%       delta   : vector containing prediction error at each trial
%       QQ      : matrix containing choice values at each trial
%
% Aroma Dabas
% October 2022
% =========================================================================


% number of trials
T = length(a);

% initialise values
q = [0.5 0.5 0.5 0.5];      % initial value of all four stimuli
PP = nan(T,2);              % probability of selecting either of the stimuli
delta = nan(1,T);           % prediction error
QQ = nan(T, size(q, 2));    % trial wise updated value
choiceProb = nan(1,T);      % trial wise choice probability

% sort presented stimuli into [HR LR]
sSorted = sort(s,2);

% loop over all trial
for t = 1:T
    
    % store the value
    QQ(t,:) = q; 
    
    % subset for the choices presented
    q_sub = q(sSorted(t,:));

    % compute choice probabilities
    p = M3_softmaxFunction(q_sub, beta);
    
    % store choice probabilities
    PP(t,:) = p;
    
    if isnan(a(t)) || a(t) == 0
        choiceProb(t) = NaN;
        delta(t) = NaN;
        
    else
        % compute choice probability for actual choice
        choiceProb(t) = p(sSorted(t,:) == a(t)); %p(a(t));

        % value update
        [q(a(t)), delta(t)] = M3_valueUpdate(alpha, q(a(t)), r(t));
    end
end

% for last trial
QQ(t+1,:) = q;

% update Q for missed trials
QQ(find(a == 0)+1,:) = [];
PP(a == 0,:) = NaN;

% compute negative log-likelihood, exluding the missed trials
NegLL = -sum(log(choiceProb(~isnan(choiceProb))));

end
