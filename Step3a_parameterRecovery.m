% In part a of the parameter recovery steps, I focus on the parameters
% recovered for RW+softmax function. The following script tests if the
% parameter used for simulating the data (here, alpha and beta parameters)
% can be recovered: Are the simulated and fit parameter values strongly
% correlated? Is there a correlation between the fitted parameters?
%
% Moreso, we also check if there is an overall bias in the estimated
% parameters, and if the variance of the estimated parameters is widely
% spread.
% 
% This code is a modification of the code provided by Bob Wilson and Anne
% Collins [2019].
%
% Modified by Aroma Dabas [dabas@cbs.mpg.de]
% 01-2020       version 1

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

% reward conditions
cond        = {'High Reward', 'Low Reward'};
n.cond      = length(cond);

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath('./HelperFunctions')
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

% color settings for plotting
AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZriver = [7, 104, 115]/256;

savePlots   = 1;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

% labels and titles
names = {'learning rate' 'softmax temperature'};
symbols = {'\alpha' '\beta'};

%% Section 2: simulate and fit data

% simulate and fit data for nRep repetitions
for count = 1:nRep

    % random parameter value for generating data
    alpha = rand;
    beta = exprnd(10);
    
    % simulate data
    [choice, reward, pt, q] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds);
    
    % fit the data
    [xf, LL] = fit_M3RescorlaWagner_v1(choice, reward, pt);

    % store simulated and fitted values, and also the - log likelihood
    fminX.sim(1,count) = alpha;
    fminX.sim(2,count) = beta;
    fminX.fit(1,count) = xf(1);
    fminX.fit(2,count) = xf(2);
    fminX.negLL(count) = LL;

    % clear repeating variables from the workspace
    clear X0 xf NegLL LB UB

end

%% Section 3: basic parameter recovey plots (simulated vs recovered)

% initiate figure
fh = figure('Name', 'Parameter Recovery'); clf;
set(gcf, 'Position', [300   513   900   370])
[~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);

% plot simulate versus fit parameters
for i= 1:size(fminX.sim,1)
    axes(ax(i)); hold on;
    plot(fminX.sim(i,:), fminX.fit(i,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
    
end

% find 'bad' alpha values
thresh = 0.25;
ind = abs(fminX.sim(1,:) - fminX.fit(1,:)) > thresh;

% mark the 'bad' alpha values
for i = 1:2
    axes(ax(i));
    plot(fminX.sim(i,ind), fminX.fit(i,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
    'markerfacecolor', [1 1 1]*0.5)
end

% set softmax parameter scale to a log scale
set(ax(1,2),'xscale', 'log', 'yscale' ,'log')

% set titles and axis labels
for i = 1:size(names, 2)
    axes(ax(i)); 
    t(i) = title(names{i});
    xlabel(sprintf('simulated %s', symbols{i}))
    ylabel(sprintf('fit %s', symbols{i})); 
end

% set font and tick sizes
set(ax, 'tickdir', 'out', 'fontsize', 18);
set(t, 'fontweight', 'normal');

% add A and B to the plots
addABCs(ax(1), [-0.07 0.08], 32);
addABCs(ax(2), [-0.1 0.08], 32, 'B');

% add reference line
set(ax, 'tickdir', 'out')
for i= 1:size(fminX.sim,1)
    axes(ax(i));
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--');
end

% save figure
if savePlots
    fh.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'Parameter_recovery', 'fit-vs-simulated.png');
    saveas(gcf, filename)
end

%% Section 4: Check if the recovered parameters are correlated

fh2 = figure('Name', 'Recovered parameters'); clf;
set(gcf, 'Position', [300   313   500   370])

% plot recovered alphas and betas
plot(fminX.fit(1,:), fminX.fit(2,:), 'o', 'color', AZriver, 'markersize', 6, 'linewidth', 1)

% set labels
xlabel(sprintf('fit %s', symbols{1}))
ylabel(sprintf('fit %s', symbols{2}));

% set font and tick sizes
set(gca, 'tickdir', 'out', 'fontsize', 18);

% save figure
if savePlots
    fh2.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'Parameter_recovery', 'fit-vs-fit.png');
    saveas(gcf, filename)
end

%% Section 3: Check for bias-variance tradeoff in parameter estimation (RMSE)

% loop over all symbols
for i = 1:length(symbols)
    
    % calculate error
    fminX.error(i,:) = fminX.sim(i,:) - fminX.fit(i,:);
    
    % calculate squared error
    fminX.sqerror(i,:) = fminX.error(i,:).^2;
    
    % take the mean of the squared error
    fminX.meansqerror(i,:) = mean(fminX.sqerror(i,:));
    
    % root mean squared error
    fminX.RMSE(i,:) = sqrt(fminX.meansqerror(i,:)); 
    
end

% print on window
fprintf('\nRMSE\nalpha = %.3f\n', fminX.RMSE(1,:))
fprintf('beta = %.3f\n', fminX.RMSE(2,:))
% --

% PLOT

% open plot and settings
figure(3); clf; hold on
set(gcf, 'Position', [500   113   400   370])       % position figure window
set(gca, 'tickdir', 'out', 'fontsize', 18);         % set tick and font size

% plot
boxplot([fminX.error(1,:)',fminX.error(2,:)'])
hold on;
line([0, 3], [0, 0], 'LineStyle', ':', 'Color', 'black');

% labels
set(gca,...
    'xticklabel', {char(945), char(946)}, 'fontsize', 14)
ylabel(sprintf('simulated - estimated\nparameter error'))
xlabel('parameters')

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Parameter_recovery', 'biasVariance_RW.png');
    saveas(gcf, filename)
end

% done