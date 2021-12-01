function [a, r, pt, PP, QQ, delta] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds, Npt)
%SIMULATE_M3RESCORLAWAGNER_V1   Function for Model 3: Rescorla Wagner with softmax function.
%
% INPUT
%       T       : total number of trials
%       alpha   : parameter alpha value
%       beta    : parameter beta value
%       rprob   : reward probability 0-1
%       rbounds : specify reward bounds as a 1X2 vector
%       Npt     : number of partial trials
%
% OUTPUT
%       a       : choices made at each trial
%       r       : reward given for each of the choices. Reward ranges from
%                 0 (unpleasant) to 1 (pleasant) with 0.5 (neutral).
%       pt      : a 1XNpt vector contatining partial trial numbers
%       PP      : choice probabilities at each trials
%       QQ      : choice values at each trial
%       delta   : prediction error at each trial
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% make list of partial trials
pt = sort(randperm(T, Npt));

% initialize variables
q  = [0.5 0.5]; % initial expected reward values
a = nan(T, 1); % action
r = nan(T, 1); % reward
PP = nan(T, size(q, 2)); % probabilities
QQ = nan(T, size(q, 2)); % value update
delta = nan(T, 1); % prediction error

% cycle over trials
for t = 1:T
    
    % store value
    QQ(t,:) = q;
    
    % compute choice probabilities
    p = M3_softmaxFunction(q, beta);
    
    % store choice probability
    PP(t,:) = p;
    
    % generate choice according to choice probabability of a_t = 2
    a(t) = choose(p(2));
    
    % determine if the choice a(t) results in a pleasant or unpleasant reward
    select = choiceReward(rprob, a(t));
%     select = binornd(1, rprob(a(t)))+1;
    
    % determine the corresponding reward value
    rpos = rewardValues(rbounds);
%     rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())];
    r(t) = rpos(select);

    % value update
    [q(a(t)), delta(t)] = M3_valueUpdate(alpha, q(a(t)), r(t), t, pt);
    
end

end
