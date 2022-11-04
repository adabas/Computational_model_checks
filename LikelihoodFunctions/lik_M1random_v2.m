function NegLL = lik_M1random_v2(a, b, s)

% LIK_M1RANDOM_v2
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a : choices vector
%       b : free parameter for estimating overall bias for one option over
%               the other
%       s   : trial wise stimuli presentation
% OUPUT:
%       NegLL : the negative log likelihood value
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

sSorted = sort(s,2);

T = length(a);

% loop over all trial
for t = 1:T
    
    % compute choice probabilities
    p = M1_randomRespondingBias(b);
    
    % compute choice probability for actual choice
    choiceProb(t) = p(sSorted(t,:) == a(t));
    %choiceProb(t) = p(a(t));
    
end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));

end
