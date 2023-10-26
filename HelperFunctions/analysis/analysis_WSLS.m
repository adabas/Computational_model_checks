function out = analysis_WSLS(a, r)
%ANALYSIS_WSLS
%
% Modified from the online script provided by Bob Wilson & Anne Collins in
% their 2019 paper "Ten Simple Rules for the computational modeling of
% behavioral data".
%
% For the LSim study, we check the probability that the agent stayed
% with a choice when the previous event for the choice was pleasant,
% neutral and unpleasant.
%
% INPUT:
%       a: 1-by-nt (number of trials) vector of choices 
%       r: 1-by-nt vector of rewards
%
% OUTPUT:
%       out: 1-by-3 vector of mean probability of staying with a choice for the
%            three event categories.
%
% Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Insitute from Human Cognitive and Brain Sciences
%
% =========================================================================
% create vector to store the previous choices
aLast = [nan a(1:end-1)];

% compare whether the previous choice is same as the current choice
stay = aLast == a;

% ---
% now split the rewards into three categories
rCategories = nan(1,length(r));
rCategories(r>=0.5) = 1;        % for high pleasantness rating
rCategories(r<0.5) = 2; %3;        % for low pleasantness rating

% create vector to store the rewards on the previous trials
rLast = [nan rCategories(1:end-1)];

% ---
% determine the mean times the participant stayed with the choice
% when the previous reward for the choice was (a) high pleasant, or (b)
%  unpleasant
pleasantStay    = nanmean(stay(rLast == 1));
unpleasantStay  = nanmean(stay(rLast == 2));

out = [unpleasantStay pleasantStay];

