% The script simulates a data set and estimates the parameters using both
% the grid search approach and the fmincon optimisers. This is a sanity
% check to see the difference between the two approaches.
%
% For a model with two free parameters, the grid search approach should be
% as good as the fmincon procedure for estimating the parameters. The
% difference should, however, arise when we go for models with more than
% two free parameters.
%
% Currently, it checks for one simulated data. In future, this can be
% implemented for a number of simulated data to get a more generalisable
% effect.
% 
% Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 
%   01-2020         version 1
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% set seed
rng(234, 'twister');

% ================== Modify ===============================================
% experiment parameters
T       = 132;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nRep    = 40;   % 100    % number of simulation repetitions
rprob   = [0.7 0.4]; % reward probability for [HR LR] stimuli
Npt     = 32;       % number of partial trials

% store information about the parameter
paramlabel  = {'alpha', 'beta'};
n.param     = length(paramlabel);

% state the bounds and bins of each parameter
bounds  = [0 1; % alpha
    0 50]; % beta
n.bin   = [20 30];

% parameter setting for simulating the data for RW+softmax function
realalpha   = .2;
realbeta    = 10;

% reward conditions
% these conditions don't seem to be used
% cond        = {'High Reward', 'Low Reward'};
% n.cond      = length(cond);

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath('./HelperFunctions')
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots   = 1;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

%% Section 2: Simulate one data set

[choice, reward, pt] = simulate_M3RescorlaWagner_v1(T, realalpha, realbeta,...
                            rprob, rbounds, Npt);

%% Section 3: fmincon optimiser

% set fmincon settings
options = optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% run optimization over 10 starting points
for count = 1:10

    % random starting points [alpha beta]
    X0  = [rand exprnd(10)];

    % create RW+softmax function handle to input in fmincon
    % You can reduce the probability for an error in the code if you avoid
    % using different functions for the simulation (simulate_M3RescorlaWagner_v1) 
    % and the likelihood (lik_M3RescorlaWagner_v1) (at least the way you
    % currently do it): The functions are very similar and share 99
    % percent of the computations. That is, if you  change anything in one
    % function you also have to change exactly the
    % same thing in the other function. This increases the probability for
    % errors, for example, if you forget to change sth. in both functions. 
    % I would either combine both function in one function with an option
    % to simulate and to compute the likelihood or to write different
    % sub-functions (e.g., softmax, delta-rule) and call them in both
    % functions to make sure both functions share exactly the same
    % sub-functions for all shared computations. 
    obFunc = @(x) lik_M3RescorlaWagner_v1(choice, reward, x(1), x(2), pt);
%   obFunc = @(x) lik_M3RescorlaWagner_v1(data.choice', data.rate.norm, X0(1), X0(2), idPartial);

    % store the lower and upper bounds of [alpha beta] parameters
    LB = [bounds(1,1) bounds(2,1)];
    UB = [bounds(1,2) bounds(2,2)];

    % run fmincon to check the best fitting parameter to the data
    [xf, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

    % store starting value (Xint), fitted value (Xfit) and the negative
    % LL values (negLL)
    fminX.initiate(1,count) = X0(1);
    fminX.initiate(2,count) = X0(2);
    fminX.fit(1,count)      = xf(1);
    fminX.fit(2,count)      = xf(2);
    fminX.negLL(count)      = NegLL;

    % clear repeating variables from the workspace
    clear X0 xf NegLL LB UB

end

% find global best
[mf,i]=min(fminX.negLL(:));
fminX.pars = fminX.fit(:,i);

%% Section 4: grid search

% create bins for both the parameters
for iParam = 1:n.param
    range = linspace(bounds(iParam,1),bounds(iParam,2),n.bin(iParam)+1);
    p{iParam} = range(2:end); % stay just off the zero bounds
end

% select each bin and determine the likelihood of the data given the
% parameters in that bin
params = nan(1,n.param);
for t = 1:(n.bin(1))
    params(1) = p{1}(t);    % alpha value
    for tt = 1:(n.bin(2))
        params(2) = p{2}(tt);   % beta value

        % determine the log likelihood value for the parameter
        llh(t,tt) = -lik_M3RescorlaWagner_v1(choice,...
            reward, params(1), params(2), pt);
    end
end

%% Section 4: Grid plot
% Plot parameters used for simulation with parameters estimated by grid search and fmincon optimiser

% create figure
fh = figure;

% find minimum and maximum values from the grid search
mi      = min(llh(:));
[ma,i]  = max(llh(:));

% create matrices containing bins of alpha and beta parameters
x = repmat(1:length(p{1,1}), length(p{1,2}), 1)';   % repeating 
y = repmat(1:length(p{1,2}), length(p{1,1}), 1);

% plot grid surface
imagesc(p{1,1}(1:end), p{1,2}(1:end), llh',[mi,ma])
colorbar

% plot the parameters
hold on
plot(realalpha, realbeta, 'xr')     % simulation parameters
plot(p{1,1}(x(i)), p{1,2}(y(i)), 'ok')  % grid search
plot(fminX.pars(1), fminX.pars(2),'*k')  % fmincon optimisation

% add labels and edit settings
xlabel('alpha')
ylabel('beta')
set(gca,'fontsize', 14)
dim = [0.65 0 0.5 0.05];
annotation('textbox',dim, 'String', '*: fmincon; o: grid search',...
    'EdgeColor', 'none', 'FontSize', 11);

% save plot
if savePlots
    fh.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Grid-vs-fmincon/Grid.png', plotFolder))
end

% done