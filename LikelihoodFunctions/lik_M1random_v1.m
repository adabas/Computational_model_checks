function NegLL = lik_M1random_v1(a, b)

% LIK_M1RANDOM_v1
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a : choices vector
%       r : reward received
%       b : free parameter for estimating overall bias for one option over
%               the other
% OUPUT:
%       NegLL : the negative log likelihood value
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   December 2019
% =========================================================================

T = length(a);

% loop over all trial
for t = 1:T
    
    % compute choice probabilities
    p = [b 1-b];
    
    % compute choice probability for actual choice
    choiceProb(t) = p(a(t));
    
end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));