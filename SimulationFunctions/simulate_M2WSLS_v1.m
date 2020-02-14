function [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, rprob)

% Function for Model 2: Noisy Win Stay Lose Sift.
%
% At a trial, agent sticks with the option that previously resulted in a
% reward. In the noisy version, the participant applies the WSLS rule with
% a probability of 1-epsilon.
%
% Reward value ranges from 0 (negative reward) to 1 (positive reward) with
% 0.5 as neutral reward.
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

% randomise 32 trials over 132 trials as partial trials
pt = sort(randperm(132,32));

for t = 1:T
    
    % compute choice probabilities
    if isnan(rLast)
        
        % first trial choose randomly
        p = [0.5 0.5];
        
    else
        
        % choice depends on last reward
        if rLast >= 0.5
            
            % win stay (with probability 1-epsilon)
            p = epsilon/2*[1 1];
            p(aLast) = 1 - epsilon / 2;
            
        else
            
            % lose shift (with probability 1-epsilon)
            p = (1-epsilon/2) * [1 1];
            p(aLast) = epsilon / 2;
            
        end
    end
    
    % make choice according to choice probababilities
    a(t) = choose(p);
    
    % determine reward as per the reward probabilities
    % select the reward probability of the selected choice
    switch a(t)
        case 1
            rprob_tmp = rprob(1);
           
        case 2
            rprob_tmp = rprob(2);
    end
    
    % determine if this choice should receive a reward or not based on the
    % selected reward probability
    x = rand;
    if x < rprob_tmp
        select = 1;  % reward
    else
        select = 2;  % no reward
    end
    
    % determine the corresponding reward value (range 0 to 1 with 0.5 as neutral)
    rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())]; 
    r(t) = rpos(select);

    aLast = a(t);
    rLast = r(t);
end
