function NegLL = lik_M2WSLS_v1(a, r, epsilon)

% LIK_M2WSLS_v1
% Function to compute the negative log-likelihood values for fitting the model to the data.
%
% INPUT:
%       a       : choices vector
%       r       : reward received
%       epsilon : probability with which to stick to rewarded stimuli
%
% OUPUT:
%       NegLL : the negative log likelihood value
%
% Aroma Dabas [dabas@cbs.mpg.de]
%   January 2020
% =========================================================================


% last reward/action (initialize as nan)
rLast = nan;
aLast = nan;

T = length(a);

% loop over all trial
for t = 1:T
    
     % compute choice probabilities
    if isnan(rLast)
        
        % first trial choose randomly
        p = [0.5 0.5];
        
    else
        
        % choice depends on last reward
        if rLast >= 0.5
            
            % win stay (with probability epsilon)
            p = epsilon/2*[1 1];
            p(aLast) = 1-epsilon/2;
            
        else
            
            % lose shift (with probability 1-epsilon)
            p = (1-epsilon/2) * [1 1];
            p(aLast) = epsilon / 2;
            
        end
    end
    
    % compute choice probability for actual choice
    choiceProb(t) = p(a(t));
    
    aLast = a(t);
    rLast = r(t);
end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));