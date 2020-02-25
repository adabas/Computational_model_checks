function [NegLL, PP, delta] = lik_M3RescorlaWagner_v1(a, r, alpha, beta, pt)

% store the initial (liking) value of the stimuli
Q = [0.5 0.5];

% number of trials
T = length(a);

% store the evolving probabilities and prediction error
PP = nan(T, 2);
delta = nan(1, T);

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
        delta(t) = 0;
    else
        delta(t)   = r(t) - Q(a(t));
    end
    Q(a(t)) = Q(a(t)) + alpha * delta(t);

end

% compute negative log-likelihood
NegLL = -sum(log(choiceProb));
