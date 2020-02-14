function [a, r, pt] = simulate_M1random_v1(T, rbounds, b, rprob)

% Function for Model 1: Random responding
%
% Model captures participants' behaviour to select random choices, perhaps
% with some overall bias for one option (HR or LR stimuli) over the other.
% 
% Reward value ranges from 0 (negative reward) to 1 (positive reward) with
% 0.5 as neutral reward.
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
%   24-10-2019      Included partial trials.
%   08-01-2020      Changed binary reward (0 and 1) to continuous reward
%                   options.
%   28-01-2020      Update rewards as per the uncorrelated reward probabilities
% =========================================================================

% randomise 32 trials over 132 trials as partial trials
pt = sort(randperm(132,32));

for t = 1:T
    
    % compute choice probability; here set to 0.5 always
    p = [b 1-b];

    % make choice according to the computed choice probability
    a(t) = choose(p);
    
    % select the reward probability of the selected choice
    switch a(t)
        case 1
            rprob_tmp = rprob(1);
           
        case 2
             rprob_tmp = rprob(2);
    end
    
    % determine if this choice should receive a reward or not based on the
    % reward probability
    x = rand;
    if x < rprob_tmp
        select = 1;  % reward
    else
        select = 2;  % no receive reward
    end
    
    % if the trial is a partial trial
    if ismember(t, pt)
        
        % find the previous trial when a(t) was selected
        tp = find(a(1:(t-1)) == a(t), 1, 'last' );
        if isempty(tp)
            r(t) = 0.5;     % neutral reward
        else
            r(t) = r(tp);   % select reward received on the previous trial for the stimuli
        end
    
    % if not a partial trial, select the degree of reward
    else
        rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())]; 
        r(t) = rpos(select);
        
    end
end

end
