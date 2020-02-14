function [a, r, pt, PP] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds)

% Function for Model 3: Rescorla Wagner with softmax function.
%
% INPUT
%       T       : total number of trials
%       alpha   : parameter alpha value
%       beta    : parameter beta value
%       rprob   : reward probability 0-1
%
% OUTPUT
%       a       : choices made at each trial
%       r       : reward given for each of the choices. Reward ranges from
%                 0 (unpleasant) to 1 (pleasant) with 0.5 (neutral).
%       q       : the value of the stimuli at the starting of the trial and
%                 at the end of all trials.
%       pt      : trial numbers that were partial, i.e. no reward was
%                 presented at this trial
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
%   24-10-2019      Included partial trials.
%   08-01-2020      Changed binary reward (0 and 1) to continuous reward
%                   options.
%   09-01-2020      Save the value of the stimuli
%                   Save partial trials (pt)
%   10-01-2020      Cleaned the input and output arguments
%   21-01-2020      Included reward probability
%   23-01-2020      Set initial value of the stimuli to 0.5
%   29-01-2020      Updated rewards as per the uncorrelated reward
%                   probabilities
% =========================================================================

% initiate the starting values
q  = [0.5 0.5];

% create storage for updating probabilities
PP = nan(T, size(q, 2));

% make list of partial trials
pt = sort(randperm(132,32));

for t = 1:T
    
    % compute choice probabilities
    ev  = exp(beta*q);
    sev = sum(ev);
    p = ev / sev;
    
    % store the trial probability
    PP(t,:) = p;
    
    % make choice according to choice probababilities
    a(t) = choose(p);
    
    % select the reward probability of the selected choice
    switch a(t)
        case 1
            rprob_tmp = rprob(1);
           
        case 2
             rprob_tmp = rprob(2);
    end
    
    % determine if this choice should receive a reward or not based on the
    % selected reward probability
    x = rand;
    if x < rprob_tmp
        select = 1;  % reward
    else
        select = 2;  % no receive reward
    end
    
    % determine the corresponding reward value
    rpos = [abs(rbounds(2)-0.5*rand()) abs(rbounds(1)-0.5*rand())]; 
    r(t) = rpos(select);
    
    % if the trial is a partial trial
    if ismember(t, pt)
        delta = 0;  % no update
    else
        delta = r(t) - q(a(t));     
    end
    
    % update the value
    q(a(t)) = q(a(t)) + alpha * delta;

end


end
