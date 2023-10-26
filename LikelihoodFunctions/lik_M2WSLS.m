function [NegLL, choiceProb, PP] = lik_M2WSLS(a, r, epsilon, s)

% LIK_M2WSLS
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       epsilon : probability with which to stick to rewarded stimuli
%       s       : trial wise stimuli presentation
%
% OUPUT:
%       NegLL   : the negative log likelihood value
%       choiceProb: estimated choice
%       PP      : trial wise choice probabilities for both HR and LR stim
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% combination options
cName = {'one', 'two', 'three', 'four'};

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

T = length(a);

% track trial wise probability for both HR and LR stimuli
PP = nan(T,2);

% loop over all trial
for t = 1:T
    
    % select last choice and reward for the trial combination
    trial_comb = s_comb{t};
    aLast = c.(sprintf('%s', trial_comb)).aLast;
    rLast = c.(sprintf('%s', trial_comb)).rLast;
    trial_choices = sSorted(t, :);
    
    % compute choice probabilities
    p = M2_WSLSprob(aLast, rLast, epsilon);
    PP(t,:) = p;
    
    % compute choice probability for the selected choice
    if a(t) == 0
        choiceProb(t) = NaN;
    else
        idChoice = find(trial_choices == a(t));
        choiceProb(t) = p(idChoice);
    
        % update last choice
        [c.(sprintf('%s', trial_comb)).aLast, c.(sprintf('%s', trial_comb)).rLast] = M2_updateChoice(idChoice, r(t));
    end

end

% update PP for missed trials
PP(a == 0,:) = NaN;

% compute negative log-likelihood
NegLL = -sum(log(choiceProb(~isnan(choiceProb))));
end