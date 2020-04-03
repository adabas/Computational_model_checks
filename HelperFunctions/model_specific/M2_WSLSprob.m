function p = M2_WSLSprob(aLast, rLast, epsilon)

% M2_WSLSPROB
% Function to compute choice probabilities based on the noisy win stay lose
% shift (noisy WSLS) model.
% 
% INPUT
%       aLast   : last choice
%       rLast   : reward for the last choice
%       epsilon : probability for sticking to the rewarded stimuli
% 
% OUTPUT
%       p       : vector containing choice probabilities
%
% Aroma Dabas [dabas@cbs.mpg.de]
% March 2020
% =========================================================================

% if this is the first trial...
if isnan(rLast)

    % ...choose randomly
    p = [0.5 0.5];

% otherwise...
else

    % ...choice depends on last reward
    if rLast >= 0.5

        % win stay (with probability epsilon)
        p = epsilon/2*[1 1];
        p(aLast) = 1 - epsilon / 2;

    else

        % lose shift (with probability 1-epsilon)
        p = (1-epsilon/2) * [1 1];
        p(aLast) = epsilon / 2;

    end
end

end