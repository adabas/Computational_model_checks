function [sim, fit] = choiceEstimation(varargin)

% choiceEstimation
% The function estimates and plots the estimated parameters trial-by-trial
% choice behaviour for the models of interest.
%
%   VARARGIN
%       ntrials     : number of trials
%       pbounds     : parameter upper and lower bounds
%       rbounds     : reward bounds
%       rewProb     : reward probability structure
%       Npt         : number of partial trials
%       nRep        : number of iterations
%       model       : model number for estimating the data
%       nMod        : number of models
%
%   OUTPUT
%       sim         : simulated data for the model
%       fit         : estimated data using the best fitting parameter
%                     values
%
% Aroma Dabas [dabas@cbs.mpg.de]
% June 2020
% =========================================================================

%% Parse input arguments

% Default values
ntrials = 100;
pbounds = [0.3 0.1 0.7 4;   % LB
    0.4 0.4 0.8 6];         % UB
rbounds = [0 1];
rewProb = [0.8 0.2];
Npt     = 0;
nRep    = 100;
model   = 1;
nMod    = 3;

% parse varargin
i = 1;
while(i<=length(varargin))
    switch varargin{i}
        case 'ntrials'
            i             = i + 1;
            ntrials       = varargin{i};
            i             = i + 1;
        case 'pbounds'
            i             = i + 1;
            pbounds       = varargin{i};
            i             = i + 1;
        case 'rbounds'
            i             = i + 1;
            rbounds       = varargin{i};
            i             = i + 1;
        case 'rewProb'
            i             = i + 1;
            rewProb       = varargin{i};
            i             = i + 1;
        case 'Npt'
            i             = i + 1;
            Npt          = varargin{i};
            i             = i + 1;
        case 'nRep'
            i             = i + 1;
            nRep          = varargin{i};
            i             = i + 1;
        case 'model'
            i             = i + 1;
            model         = varargin{i};
            i             = i + 1;
        case 'nMod'
            i             = i + 1;
            nMod          = varargin{i};
            i             = i + 1;
    end
end


%% Estimate

% iterate
for count = 1:nRep
    
    % simulate data
    if model == 1
        b = pbounds(1,1) + (pbounds(2,1)-pbounds(1,1)) .* rand(1,1);   % constrained parameter
        [a, r, pt] = simulate_M1random_v1(ntrials, rbounds, b, rewProb, Npt);
    elseif model == 2
        epsilon = pbounds(1,2) + (pbounds(2,2)-pbounds(1,2)) .* rand(1,1);
        [a, r, pt] = simulate_M2WSLS_v1(ntrials, rbounds, epsilon, rewProb, Npt);
    elseif model == 3
        alpha = pbounds(1,3) + (pbounds(2,3)-pbounds(1,3)) .* rand(1,1);
        beta = pbounds(1,4) + (pbounds(2,4)-pbounds(1,4)) .* rand(1,1);
        [a, r, pt] = simulate_M3RescorlaWagner_v1(ntrials, alpha, beta, rewProb, rbounds, Npt);
    elseif model == 4
        alpha = pbounds(1,3) + (pbounds(2,3)-pbounds(1,3)) .* rand(1,1);
        beta = pbounds(1,4) + (pbounds(2,4)-pbounds(1,4)) .* rand(1,1);
        alpha_c = pbounds(1,5) + (pbounds(2,5)-pbounds(1,5)) .* rand(1,1);
        beta_c = pbounds(1,6) + (pbounds(2,6)-pbounds(1,6)) .* rand(1,1);
        [a, r, pt] = simulate_M4RWCK_v1(ntrials, alpha, beta, alpha_c, beta_c, rewProb, rbounds, Npt);
    elseif model == 5
        alpha_c = pbounds(1,5) + (pbounds(2,5)-pbounds(1,5)) .* rand(1,1);
        beta_c = pbounds(1,6) + (pbounds(2,6)-pbounds(1,6)) .* rand(1,1);
        [a, r, pt] = simulate_M5CK_v1(ntrials, alpha_c, beta_c, rewProb, rbounds, Npt);
    end

    % store simulated data
    sim(1).a(:,count) = a;  
    sim(1).r(:,count) = r;

    % estimate the best fitting parameters for all models
    [~,~,~, pars] = fit_all_v1(a, r, pt, nMod, pbounds, a);

    % use parameter values to estimate the corresponding choice behaviour
    % - Model 1
    [a, r] = simulate_M1random_v1(ntrials, rbounds, pars(1,1), rewProb, Npt);
    fit(1).a(:,count) = a;
    fit(1).r(:,count) = r;

    % - Model 2
    [a, r] = simulate_M2WSLS_v1(ntrials, rbounds, pars(2,1), rewProb, Npt);
    fit(2).a(:,count) = a;
    fit(2).r(:,count) = r;

    % - Model 3
    [a, r] = simulate_M3RescorlaWagner_v1(ntrials, pars(3,1), pars(3,2), rewProb, rbounds, Npt);
    fit(3).a(:,count) = a;
    fit(3).r(:,count) = r;
    
    % - Model 4
    [a, r] = simulate_M4RWCK_v1(ntrials, pars(4,1), pars(4,2), pars(4,3), pars(4,4), rewProb, rbounds, Npt);
    fit(4).a(:,count) = a;
    fit(4).r(:,count) = r;
    
    % - Model 5
    [a, r] = simulate_M5CK_v1(ntrials, pars(5,1), pars(5,2), rewProb, rbounds, Npt);
    fit(5).a(:,count) = a;
    fit(5).r(:,count) = r;

end

end