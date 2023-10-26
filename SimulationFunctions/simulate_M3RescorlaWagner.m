function [a, r, s, PP, QQ, delta] = simulate_M3RescorlaWagner(T, alpha, beta, rprob, rbounds)
%SIMULATE_M3RESCORLAWAGNER_V2
% Function for Model 3: Rescorla Wagner with softmax function.
%
% INPUT
%       T       : total number of trials
%       alpha   : parameter alpha value
%       beta    : parameter beta value
%       rprob   : reward probability 0-1
%       rbounds : specify reward bounds as a 1X2 vector
%
% OUTPUT
%       a       : choices made at each trial
%       r       : reward given for each of the choices. Reward ranges from
%                 0 (unpleasant) to 1 (pleasant) with 0.5 (neutral).
%       PP      : choice probabilities at each trials
%       QQ      : choice values at each trial
%       delta   : prediction error at each trial
%       s       : trial wise stimuli presentation
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% initialize variables
q  = [0.5 0.5 0.5 0.5]; % initial expected reward values
a = nan(T, 1); % action
r = nan(T, 1); % reward
PP = nan(T, 2); % probabilities
QQ = nan(T, size(q, 2)); % value update
delta = nan(T, 1); % prediction error

% load a stimuli presentation
s = stimuliPresentation(T);
sSorted = sort(s, 2);

% cycle over trials
for t = 1:T
    
    % store value
    QQ(t,:) = q;
    
    % choices presented at trial t
    q_sub = q(sSorted(t,:));
    
    % compute choice probabilities
    p = M3_softmaxFunction(q_sub, beta);
    
    % store choice probability
    PP(t,:) = p;
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = sSorted(t, choose(p(2)));
    
    % UPDATE: for fMRI analysis, we are binaring the rewards.
    r(t) = rand < rprob(a(t));

    % value update
    [q(a(t)), delta(t)] = M3_valueUpdate(alpha, q(a(t)), r(t));
    
end

end
