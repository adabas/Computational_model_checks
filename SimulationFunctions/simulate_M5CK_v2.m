function [a, r, pt, s, PP, CK] = simulate_M5CK_v2(T, alpha_c, beta_c, rprob, rbounds, Npt)

% SIMULATE_M5CK_V2
% Function simulates data for a choice kernel model.
% Here, we compute the frequency of selecting a choice.
% 
% INPUT
%       T       : total number of trials
%       alpha_c : choice kernel alpha parameter value
%       beta_c  : choice kernel beta parameter value
%       rprob   : reward probability 0-1
%       rbounds : specify reward bounds as a 1X2 vector
%       Npt     : number of partial trials
%
% OUTPUT
%       a       : choices made at each trial
%       r       : reward given for each of the choices. Reward ranges from
%                 0 (unpleasant) to 1 (pleasant) with 0.5 (neutral).
%       pt      : a 1XNpt vector contatining partial trial numbers
%       s       : trial wise stimuli presented
%       PP      : choice probabilities at each trials
%       CK      : choice kernel update at each trial
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% make list of partial trials
pt = sort(randperm(T, Npt));

% initialise variables
k = [0 0 0 0];     % initialise choice kernel
a = nan(T, 1); % action
r = nan(T, 1); % reward
PP = nan(T, 2); % probabilities
CK = nan(T, size(k, 2)); % choice kernal update

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
                
    % generate choice according to choice probabability of a_t = 2
    a(t) = sSorted(t, choose(p(2)));
    
    % UPDATE: for fMRI analysis, we are binaring the rewards.
    r(t) = rand < rprob(a(t));

%    % determine if the choice a(t) results in a pleasant or unpleasant reward
%     select = choiceReward(rprob, a(t));
%     
%     % determine the corresponding reward value
%     rpos = rewardValues(rbounds);
%     r(t) = rpos(select);
    
    % update value and choice kernel
    k(a(t)) = M5_CKUpdate(alpha_c, k(a(t)));
    
end

end
