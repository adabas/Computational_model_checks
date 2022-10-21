function [a, r, pt, s, PP, QQ, CK, delta] = simulate_M4RWCK_v2(T, alpha, beta, alpha_c, beta_c, rprob, rbounds, Npt)

% SIMULATE_M4RWCK_V2
% Function simulates data for a Rescorla Wagner and a choice kernel model.
% Here, we compute the value update while keeping a track of how frequently
% a choice is selected.
% 
% INPUT
%       T       : total number of trials
%       alpha   : alpha parameter value
%       beta    : beta parameter value
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
%       QQ      : choice value update at each trial
%       CK      : choice kernel update at each trial
%       delta   : prediction error at each trial
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% =========================================================================

% make list of partial trials
pt = sort(randperm(T, Npt));

% initialise variables
q  = [0.5 0.5 0.5 0.5]; % initial expected values
k = [0 0 0 0];     % initialise choice kernel
a = nan(T, 1); % action
r = nan(T, 1); % reward
PP = nan(T, 2); % probabilities
QQ = nan(T, size(q, 2)); % value update
CK = nan(T, size(q, 2)); % choice kernal update
delta = nan(T, 1); % prediction error

% load a stimuli presentation
s = stimuliPresentation(T);

for t = 1:T
    
    % store value and kernel
    QQ(t,:) = q;
    CK(t,:) = k;
    
    % choices presented at trial t
    q_sub = q(s(t,:));
    k_sub = k(s(t,:));
    
    % compute choice probabilities
    p = M4_softmaxRWCK(q_sub, k_sub, beta, beta_c);
    
    % store choice probability
    PP(t,:) = p;
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = s(t, choose(p(2)));
    
    % UPDATE: for fMRI analysis, we are binaring the rewards.
    r(t) = rand < rprob(a(t));
    
%     % generate choice according to choice probabability of a_t = 2
%     a(t) = choose(p(2));
%     
%     % determine if the choice a(t) results in a pleasant or unpleasant reward
%     select = choiceReward(rprob, a(t));
%     
%     % determine the corresponding reward value
%     rpos = rewardValues(rbounds);
%     r(t) = rpos(select);
    
    % update value and choice kernel
    [q(a(t)), k(a(t)), delta(t)] = M4_valueUpdate(alpha, alpha_c, q(a(t)), k(a(t)), r(t), t, pt);         
    
end

end
