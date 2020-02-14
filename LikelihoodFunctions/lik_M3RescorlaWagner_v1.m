function [NegLL, PP] = lik_M3RescorlaWagner_v1(a, r, alpha, beta, pt)

% store the initial (liking) value of the stimuli
Q = [0.5 0.5];

% number of trials
T = length(a);

% store the evolving probabilities
PP = nan(T, 2);

% loop over all trial
for t = 1:T
    
    % compute choice probabilities
    ev  = exp(beta*Q);
    sev = sum(ev);
    p   = ev / sev;
    
    % store choice probabilities
    PP(t,:) = p;
    
    % compute choice probability for actual choice
    choiceProb(t) = p(a(t));
    
    % update values but not for partial trials
    if ismember(t, pt)
        delta = 0;
    else
        delta   = r(t) - Q(a(t));
    end
    Q(a(t)) = Q(a(t)) + alpha * delta;

end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
