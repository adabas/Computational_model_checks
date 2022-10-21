function NegLL = lik_M2WSLS_v2(a, r, epsilon, s)

% LIK_M2WSLS_v2
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       epsilon : probability with which to stick to rewarded stimuli
%       s       : trial wise stimuli presentation
%
% OUPUT:
%       NegLL : the negative log likelihood value
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

% loop over all trial
for t = 1:T
    
    % select for the trial combination
    trial_comb = s_comb{t};
    aLast = c.(sprintf('%s', trial_comb)).aLast;
    rLast = c.(sprintf('%s', trial_comb)).rLast;
    trial_choices = sSorted(t, :);
    
    % compute choice probabilities
    p = M2_WSLSprob(aLast, rLast, epsilon);
    
    % compute choice probability for actual choice
    idChoice = find(sSorted(t,:) == a(t));
    choiceProb(t) = p(idChoice);%p(a(t));
    
    % update last choice
    %[aLast, rLast] = M2_updateChoice(a(t), r(t));
    [c.(sprintf('%s', trial_comb)).aLast, c.(sprintf('%s', trial_comb)).rLast] = M2_updateChoice(idChoice, r(t));

end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
end