function [a, r, pt, s] = simulate_M2WSLS_v2(T, rbounds, epsilon, rprob, Npt)
%SIMULATE_M2WSLS_V3   Function for Model 2: Noisy Win Stay Lose Sift.
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
% October 2022
% =========================================================================

% combination options
cName = {'one', 'two', 'three', 'four'};

% create trial specific stimuli presentation
s = stimuliPresentation(T);

% sort into HR vs LR
sSorted = sort(s,2);

s_comb = cell(96,1);
s_comb(sSorted(:,1) == 1 & sSorted(:,2) == 3) = cName(1);
s_comb(sSorted(:,1) == 1 & sSorted(:,2) == 4) = cName(2);
s_comb(sSorted(:,1) == 2 & sSorted(:,2) == 3) = cName(3);
s_comb(sSorted(:,1) == 2 & sSorted(:,2) == 4) = cName(4);

% last reward/action (initialize as nan) for each stimuli combination pair
for i = 1:numel(cName)
    c.(sprintf('%s', cName{i})).rLast = nan;
    c.(sprintf('%s', cName{i})).aLast = nan;
end

% randomise Npt trials over T trials as partial trials
pt = sort(randperm(T,Npt));

% initialize variables
a = nan(T, 1); % action
r = nan(T, 1); % reward

% cycle over trials
for t = 1:T
    
    % select for the trial combination
    trial_comb = s_comb{t};
    aLast = c.(sprintf('%s', trial_comb)).aLast;
    rLast = c.(sprintf('%s', trial_comb)).rLast;
    trial_choices = sSorted(t, :);
    
    % compute choice probabilities
    p = M2_WSLSprob(aLast, rLast, epsilon);
    
    % generate choice according to choice probabability of a_t = 2
    idChoice = choose(p(2));
    a(t) = trial_choices(idChoice);
    
    % UPDATE: for fMRI analysis, we are binaring the rewards.
    r(t) = rand < rprob(a(t));
    
%     % determine if the choice a(t) results in a pleasant or unpleasant reward
%     select = choiceReward(rprob, a(t));
% %     select = binornd(1, rprob(a(t)))+1;
%     
%     % determine the corresponding reward value (range 0 to 1 with 0.5 as neutral)
%     rpos = rewardValues(rbounds);
% %     rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
%     r(t) = rpos(select);
    
    % update last choice trial
    %[aLast, .rLast] = M2_updateChoice(a(t), r(t));
    [c.(sprintf('%s', trial_comb)).aLast, c.(sprintf('%s', trial_comb)).rLast] = M2_updateChoice(idChoice, r(t));

end

end
