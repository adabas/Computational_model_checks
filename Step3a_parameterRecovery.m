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
nRep    = 100; %1000;     % number of simulation repetitions
nFit    = 10;       % iterate over the fitting function for multiple starting points
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
    beta = exprnd(4);

    % simulate data
    [choice, reward, pt, q] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds, Npt);
    
    % iterate the best fit
    for iter = 1:nFit
        
        % fit the data
        [xf, LL] = fit_M3RescorlaWagner_v1(choice, reward, pt);
        
        fit(1,iter) = xf(1);
        fit(2,iter) = xf(2);
        negLL(iter) = -LL;  % store the negative log likelihood
    end
    
    % find global best
    [mf,i]=min(negLL(:));
    pars = fit(:,i);
    
    % Here you could add the option to cycle over multiple (random) starting points.
    % I.e., call [xf, LL] = fit_M3RescorlaWagner_v1(choice, reward, pt); 
    % n times and select parameter of the iteration with the best fit (lowest neg. log-likelihood)
    % At the moment you only use one random starting point. Using more
    % starting point might improve the estimates in some cases. 
    % Especially in models with more free parameters, you'll notice a
    % difference. Right that was a point we discussed. Thanks for the
    % reminder!

    % store simulated and fitted values, and also the log likelihood
    fminX.sim(1,count) = alpha;
    fminX.sim(2,count) = beta;
    fminX.fit(1,count) = pars(1);
    fminX.fit(2,count) = pars(2);
    fminX.negLL(count) = mf;

%     % clear repeating variables from the workspace
%     clear xf NegLL LL mf i 

end

%% Section 3: basic parameter recovery plots (simulated vs recovered)

% initiate figure
fh = figure('Name', 'Parameter Recovery'); clf;
set(gcf, 'Position', [300   513   900   370])
[~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);

% plot simulate versus fit parameters
for i = 1:size(fminX.sim,1)
    axes(ax(i)); hold on;
    plot(fminX.sim(i,:), fminX.fit(i,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
end

% find 'bad' parameter values
% I think in the previous version, you also used the "bad" alpha parameters
% for the plot of the "bad" beta values. I changed this and added a
% different threshold for beta. 

% mark the 'bad' parameter values
for i = 1:2
    
    if i == 1
        thresh = 0.25;
    else
        thresh = 5;
    end
    ind = abs(fminX.sim(i,:) - fminX.fit(i,:)) > thresh;
    axes(ax(i));
    plot(fminX.sim(i,ind), fminX.fit(i,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
    'markerfacecolor', [1 1 1]*0.5)
end

% set softmax parameter scale to a log scale
% Not sure if this helps me. I think it's not necessary to see if it works, especially with lower beta values.
% But of course not wrong. Yes, in this case doesn't seem necessary. But
% I'll keep this commented out.
% set(ax(1,2),'xscale', 'log', 'yscale' ,'log')

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

% plot settings: set y axis limits for the parameters
yLim = [[-1 1]; [-20 20]];

% initiate figure
fh3 = figure('Name', 'Parameter Estimation Error'); clf;
set(gcf, 'Position',[500   113   550   370])
%set(gca, 'tickdir', 'out', 'fontsize', 18);         % set tick and font size
[~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);

% plot estimation error
for i = 1:size(fminX.error,1)
    axes(ax(i)); hold on;       % select axis
    boxplot(fminX.error(i,:)', 'Labels',{[]})  % boxplot
    hold on;
    ylim(yLim(i,:))
    line([0, 3], [0, 0], 'LineStyle', ':', 'Color', 'black'); % add reference line
end

% set titles and axis labels
for i = 1:size(names, 2)
    axes(ax(i));                % select axis
%    t(i) = title(names{i});
    xlabel(sprintf('%s', symbols{i}))
    ylabel(sprintf('simulated - estimated\nparameter error'));
end

% set font and tick sizes
set(ax, 'tickdir', 'out', 'fontsize', 12);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Parameter_recovery', 'biasVariance_RW.png');
    saveas(gcf, filename)
end

% done