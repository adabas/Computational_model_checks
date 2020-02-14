% Modified from the online script provided by Robert Wilson & Anne Collins in
% their 2019 paper "Ten Simple Rules for the computational modeling of
% behavioral data".
%
% For the LSim study, we check the probability that the agent stayed
% with a choice when the previous event for the choice was pleasant,
% neutral and unpleasant.
%
% INPUT:
%       1-by-nt (number of trials) vector of choices (a)
%       1-by-nt vector of rewards (r)
%
% OUTPUT:
%       1-by-3 vector of mean probability of staying with a choice for the
%       three event categories.
%
% Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Insitute from Human Cognitive and Brain Sciences
%
%   version 1:  12-2019     Split rewards into the three event categories.
%
% =========================================================================

function out = analysis_WSLS_v1(a, r)

% create vector to store the previous choices
aLast = [nan a(1:end-1)];

% compare whether the previous choice is same as the current choice
stay = aLast == a;

% ---
% now split the rewards into three categories
rCategories = nan(1,length(r));
rCategories(r>=0.6) = 1;        % for high pleasantness rating
rCategories(r>0.4 & r<0.6) = 2; % for neutral pleasantness rating
rCategories(r<=0.4) = 3;        % for low pleasantness rating

% create vector to store the rewards on the previous trials
rLast = [nan rCategories(1:end-1)];

% ---
% determine the mean times the participant stayed with the choice
% when the previous reward for the choice was (a) high pleasant, (b)
% neutral or (c) unpleasant
pleasantStay    = nanmean(stay(rLast == 1));
neutralStay     = nanmean(stay(rLast == 2));
unpleasantStay  = nanmean(stay(rLast == 3));

out = [unpleasantStay neutralStay pleasantStay];

