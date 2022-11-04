function [a, r, pt, s] = simulate_M1random_v2(T, rbounds, b, rprob, Npt)
%SIMULATE_M1RANDOM_V2   Function for Model 1: Random responding
%
% Model captures participants' behaviour to select random choices, perhaps
% with some overall bias for one option (left or right button press) over the other.
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
%           s       :   trial wise stimuli presentation
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% create trial specific stimuli presentation
s = stimuliPresentation(T); 
sSorted = sort(s,2); % [HR LR] sorted

% randomise Npt trials over T trials as partial trials
pt = sort(randperm(T, Npt));

% initialize variables
a = nan(T, 1); % action
r = nan(T, 1); % reward

% cycle over trials
for t = 1:T
    
    % determine choice probabilites based on bias parameter b
    p = M1_randomRespondingBias(b);
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = sSorted(t, choose(p(2)));
    
    % reward
    r(t) = rand < rprob(a(t));
    
%     % determine if the choice a(t) results in a pleasant or unpleasant reward
%     select = choiceReward(rprob, a(t));
% %     select = binornd(1, rprob(a(t)))+1;
%     
%     % if the trial is a partial trial...
%     if ismember(t, pt)
%         
%         % ...find the previous trial when a(t) was selected
%         tp = find(a(1:(t-1)) == a(t), 1, 'last');
%         if isempty(tp)
%             r(t) = 0.5;     % neutral reward
%         else
%             r(t) = r(tp);   % select reward received on the previous trial for the stimulus
%         end
%         
%     else
%         rpos = rewardValues(rbounds);
% %        rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
%         r(t) = rpos(select);
%         
%     end
end

end
