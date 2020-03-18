function [a, r, pt] = simulate_M1random_v1(T, rbounds, b, rprob, Npt)
%SIMULATE_M1RANDOM_V1   Function for Model 1: Random responding
%
% Model captures participants' behaviour to select random choices, perhaps
% with some overall bias for one option (HR or LR stimuli) over the other.
%
% Reward value ranges from 0 (unpleasant reward) to 1 (pleasant reward) with
% 0.5 as neutral reward.
%
% INPUT :
%           T       :   total number of trials
%           rbounds :   specify reward bounds as a 1X2 vector
%           b       :   option bias parameter value ranging from 0 to 1
%                       with 0.5 indicating no bias.
%           rprob   :   reward probabilities for each stimuli as a 1X2 vector
%           Npt     :   number of partial trials
%
% OUTPUT :
%           a       :   TX1 vector indicating choices at each trial
%           r       :   TX1 vector indicating rewards received for choices
%                       at each trial
%           pt      :   a 1XNpt vector containing partial trial numbers
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
%   24-10-2019      Included partial trials.
%   08-01-2020      Changed binary reward (0 and 1) to continuous reward
%                   options.
%   28-01-2020      Update rewards as per the uncorrelated reward probabilities
% =========================================================================

% randomise Npt trials over T trials as partial trials
pt = sort(randperm(T, Npt));

% initialize variables
a = nan(T, 1); % action
r = nan(T, 1); % reward

% cycle over trials
for t = 1:T
    
    % determine choice probabilites based on bias parameter b
    p = [b 1-b];
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = choose(p(2));
    
    % determine if the choice a(t) results in a pleasant or unpleasant reward
    select = binornd(1, rprob(a(t)))+1;
    
    % if the trial is a partial trial...
    if ismember(t, pt)
        
        % ...find the previous trial when a(t) was selected
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
        % Answer: Actually not because of the partial trials. But because
        % the task rewards are not binary but varies on a continuous scale.
        % We could split the scale to make the reward binary but I wonder
        % if that is over simplification. Ok, we'll keep this in mind and
        % discuss after I've reviewed the rest.
        
    else
        rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
        r(t) = rpos(select);
        
    end
end

end
