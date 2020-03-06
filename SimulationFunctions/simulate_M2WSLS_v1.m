function [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, rprob, Npt)
%SIMULATE_M2WSLS_V1   Function for Model 2: Noisy Win Stay Lose Sift.
%
% At a trial, agent sticks with the option that previously resulted in a
% reward. In the noisy version, the participant applies the WSLS rule with
% a probability of 1-epsilon.
%
% Reward value ranges from 0 (unpleasant reward) to 1 (pleasant reward) with
% 0.5 as neutral reward.
%
% INPUT:
%       T       :   total number of trials
%       rbounds :   specify reward bounds as a 1X2 vector to generate reward values
%       epsilon :   probability with which to stick to rewarded stimuli
%       rprob   :   reward probabilities for each stimuli as a 1X2 vector
%       Npt     :   number of partial trials
%
% OUTPUT:
%       a       :   TX1 vector indicating choices at each trial
%       r       :   TX1 vector indicating rewards received for choices
%                   at each trial
%       pt      :   a 1XNpt vector containing partial trial numbers
%
%
% Originally written by Rob Wilson & Anne Collins (2018).
%
% Modified by Aroma Dabas.
% Version 2     29-01-2020      Included probabilistic rewards
%
% =========================================================================

% last reward/action (initialize as nan)
rLast = nan;
aLast = nan;

% randomise 32 trials over T trials as partial trials
pt = sort(randperm(T,Npt));

% initialize variables
a = nan(T, 1); % action
r = nan(T, 1); % reward

% cycle over trials
for t = 1:T
    
    % compute choice probabilities
    if isnan(rLast)
        
        % first trial choose randomly
        p = [0.5 0.5];
        
    else
        
        % choice depends on last reward
        if rLast >= 0.5
            
            % win stay (with probability epsilon)
            p = epsilon/2*[1 1];
            p(aLast) = 1 - epsilon / 2;
            
        else
            
            % lose shift (with probability 1-epsilon)
            p = (1-epsilon/2) * [1 1];
            p(aLast) = epsilon / 2;
            
        end
    end
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = choose(p(2));
    
    % determine if the choice a(t) results in a pleasant or unpleasant reward
    select = binornd(1, rprob(a(t)))+1;
    
    % determine the corresponding reward value (range 0 to 1 with 0.5 as neutral)
    rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
    r(t) = rpos(select);
    
    aLast = a(t);
    rLast = r(t);
end

end
