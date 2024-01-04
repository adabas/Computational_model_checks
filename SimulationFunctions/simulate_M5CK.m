function [a, r, s, PP, CK] = simulate_M5CK(T, alpha_c, beta_c, rprob, rbounds)

% SIMULATE_M5CK
% Function simulates data for a choice kernel model.
% Here, we compute the frequency of selecting a choice.
% 
% INPUT
%       T       : total number of trials
%       alpha_c : choice kernel alpha parameter value
%       beta_c  : choice kernel beta parameter value
%       rprob   : reward probability 0-1
%       rbounds : specify reward bounds as a 1X2 vector
%
% OUTPUT
%       a       : choices made at each trial
%       r       : reward given for each of the choices. Reward ranges from
%                 0 (unpleasant) to 1 (pleasant) with 0.5 (neutral).
%       s       : trial wise stimuli presented
%       PP      : choice probabilities at each trials
%       CK      : choice kernel update at each trial
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% initialise variables
k = [0 0 0 0];     % initialise choice kernel
a = nan(T, 1); % action
r = nan(T, 1); % reward
PP = nan(T, 2); % probabilities
CK = nan(T, size(k, 2)); % choice kernel update

% load a stimuli presentation
s = stimuliPresentation(T);
sSorted = sort(s,2);

for t = 1:T
    
    % store value
    CK(t,:) = k;
    
    % choices presented at trial t
    k_sub = k(sSorted(t,:));
    
    % compute choice probabilities
    p = M5_softmaxCK(k_sub, beta_c);
    
    % store choice probability
    PP(t,:) = p;
                
    % generate choice according to choice probability of a_t = 2
    a(t) = sSorted(t, choose(p(2)));
    
    % generate binarized reward
    r(t) = rand < rprob(a(t));
    
    % update value and choice kernel
    k(a(t)) = M5_CKUpdate(alpha_c, k(a(t)));
    
end

end
