function [a, r, pt] = simulate_M1random_v1(T, rbounds, b, rprob, Npt)
%SIMULATE_M1RANDOM_V1   Function for Model 1: Random responding
%
%   ADD INPUT AND OUTPUT.
%
% Model captures participants' behaviour to select random choices, perhaps
% with some overall bias for one option (HR or LR stimuli) over the other.
%
% Reward value ranges from 0 (unpleasant reward) to 1 (pleasant reward) with
% 0.5 as neutral reward.
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
%   24-10-2019      Included partial trials.
%   08-01-2020      Changed binary reward (0 and 1) to continuous reward
%                   options.
%   28-01-2020      Update rewards as per the uncorrelated reward probabilities
% =========================================================================

% randomise 32 trials over T trials as partial trials
pt = sort(randperm(T, Npt));

% initialize variables
a = nan(length(pt), 1); % action
r = nan(length(pt), 1); % reward

% cycle over trials
for t = 1:T
    
    % determine choice probabilites based on bias parameter b
    p = [b 1-b];
    
    % gererate choice according to choice probabability of a_t = 2
    a(t) = choose(p(2));
    
    % select the reward probability of the selected choice
    % Is there any particular reason for this? Does my simplification below work
    % for you?
    %     switch a(t)
    %         case 1
    %             rprob_tmp = rprob(1);
    %
    %         case 2
    %              rprob_tmp = rprob(2);
    %     end
    
    % determine if this choice should receive a reward or not based on the
    % reward probability
    %     x = rand;
    %     if x < rprob_tmp
    %         select = 1;  % reward
    %     else
    %         select = 2;  % no receive reward
    %     end
    select = binornd(1, rprob(a(t)))+1;
    
    % if the trial is a partial trial
    if ismember(t, pt)
        
        % find the previous trial when a(t) was selected
        tp = find(a(1:(t-1)) == a(t), 1, 'last');
        if isempty(tp)
            r(t) = 0.5;     % neutral reward
        else
            r(t) = r(tp);   % select reward received on the previous trial for the stimulus
        end
        
        % if not a partial trial, select the degree of reward
        % Question: Why don't you just provide 0 and 1 as reward? Probably
        % because of the partial trials, but I wonder if we could simplify
        % this.
    else
        rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
        r(t) = rpos(select);
        
    end
end

end
