function [NegLL, choiceProb, PP] = lik_M1random_v2(a, b, s)

% LIK_M1RANDOM_v2
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a : choices vector
%       b : free parameter for estimating overall bias for one option over
%               the other
%       s   : trial wise stimuli presentation
% OUPUT:
%       NegLL   : the negative log likelihood value
%       choiceProb: choice probability for the selected stimuli
%       PP      : probability of selecting HR and LR stimuli
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

sSorted = sort(s,2);

T = length(a);

% track probability of selecting a HR and LR stimuli
PP = nan(T,2);

% loop over all trial
for t = 1:T
    
    % compute choice probabilities
    p = M1_randomRespondingBias(b);
    PP(t,:) = p;
    
    % compute choice probability for actual choice
    choiceProb(t) = p(sSorted(t,:) == a(t));
    %choiceProb(t) = p(a(t));
    
end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));

end
