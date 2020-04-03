% In part b of the parameter recovery steps, I focus on the parameters
% recovered for all the models of interest. The following script tests the
% bias-variance trade off: is there is an overall bias in the estimated
% parameters, and if the variance of the estimated parameters is widely
% spread.
% 
% This code is a modification of the code provided by Rob Wilson and Anne
% Collins [2019].
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
%   01-2020     Version 1
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% set seed
rng(234, 'twister');

% ================== Modify ===============================================
% experiment parameters
T       = 132;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nRep    = 1000;     % number of simulation repetitions
rprob   = [0.7 0.4]; % reward probability for [HR LR] stimuli
Npt     = 32;       % number of partial trials

% reward conditions
cond        = {'High Reward', 'Low Reward'};
n.cond      = length(cond);

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath(genpath('./HelperFunctions'))
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

% some color settings for plotting
global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

savePlots   = 1;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

% labels and titles
names = {'learning rate' 'softmax temperature', 'option bias', 'WSLS probabilistic'};
symbols = {'\alpha' '\beta' 'b' '\epsilon'};

%% Section 2a: simulate and fit data for RW + softmax function

% simulate and fit data for nRep repetitions
for count = 1:nRep

    % random values for generating data
    alpha = rand;
    beta = exprnd(10);
    
    % simulate data
    [choice, reward, pt, q] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds, Npt);
    
    % fit the data
    [xf, LL] = fit_M3RescorlaWagner_v1(choice, reward, pt);

    % store simulated values, fitted values, and - loglikelihood
    fminX.sim(1,count) = alpha;
    fminX.sim(2,count) = beta;
    fminX.fit(1,count) = xf(1);
    fminX.fit(2,count) = xf(2);
    fminX.negLL(count) = LL;

    % clear repeating variables from the workspace
    clear xf LL choice reward pt

end

%% Section 2b: simulate and fit data for random responding model


for count = 1:nRep
    
    % random parameter value for generating data
    b = rand;
    
    % simulate data
    [choice] = simulate_M1random_v1(T, rbounds, b, rprob, Npt);
    
    % fit the data
    [xf, LL] = fit_M1random_v1(choice);
    
    % store simulated values, fitted values and the - loglikelihood
    fminX.sim(3,count) = b;
    fminX.fit(3,count) = xf;
    fminX.negLL(2,count) = LL;
    
    clear xf LL choice
end


%% Section 2c: simulate and fit data for noisy WSLS model

for count = 1:nRep
    
    % random parameter value for generating data
    epsilon = rand;
    
    % simulate data
    [choice, reward] = simulate_M2WSLS_v1(T, rbounds, epsilon, rprob, Npt);
    
    % fit the data
    [xf, LL] = fit_M2WSLS_v1(choice, reward);
    
    % store simulated values, fitted values and the - loglikelihood
    fminX.sim(4,count) = b;
    fminX.fit(4,count) = xf;
    fminX.negLL(3,count) = LL;
    
    clear xf LL choice reward
end

%% Section 3: Check for bias-variance tradeoff in parameter estimation (RMSE)

% loop over the free parameters
for i = 1:length(symbols)
    
    % calculate the error between the estimated and simulated parameter
    fminX.error(i,:) = fminX.sim(i,:) - fminX.fit(i,:);
    
    % calculate squared error
    fminX.sqerror(i,:) = fminX.error(i,:).^2;
    
    % take the mean of the squared error
    fminX.meansqerror(i,:) = mean(fminX.sqerror(i,:));
    
    % root mean squared error
    fminX.RMSE(i,:) = sqrt(fminX.meansqerror(i,:));
    
end

% print on window
fprintf('\nRMSE values\nalpha = %.3f\n', fminX.RMSE(1,:))
fprintf('beta = %.3f\n', fminX.RMSE(2,:))
fprintf('b = %.3f\n', fminX.RMSE(3,:))
fprintf('epsilon = %.3f\n', fminX.RMSE(4,:))
% --

% PLOT

% open figure and set settings
figure(3); clf; hold on
set(gcf, 'Position', [500   113   600   370])       % position figure window
set(gca, 'tickdir', 'out', 'fontsize', 18);         % set tick and font size

% plot
boxplot([fminX.error(1,:)', fminX.error(2,:)', fminX.error(3,:)', fminX.error(4,:)'])
line([0, 5], [0, 0], 'LineStyle', ':', 'Color', 'black');

% labels
xlabel('parameters')
ylabel(sprintf('simulated - estimated\nparameter error'))
set(gca,...
    'xticklabel', {char(945), char(946), 'b', char(949)}, 'fontsize', 14)

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Parameter_recovery', 'biasVariance_allModels.png');
    saveas(gcf, filename)
end

% done
