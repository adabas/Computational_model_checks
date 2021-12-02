% The script tests parameter recovery for the models: simple RW, random
% responding, WSLS and RW-CK. This plots:
% 1. Scatter plot with fitted vs simulated parameters
% 2. Correlation between parameters of a) RW model and b) RW-CK model
% 3. bias-variance trade off: is there is an overall bias in the estimated
%    parameters, and if the variance of the estimated parameters is widely
%    spread.
% 
% This code is a modification of the code provided by Bob Wilson and Anne
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
T       = 100;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nRep    = 100;     % number of simulation repetitions
rprob   = [0.8 0.35]; % reward probability for [HR LR] stimuli
Npt     = 0;       % number of partial trials
%pbounds = [0 0 0 0 0 0;
%    1 1 1 7 1 7];    % bounds [lower; upper] * parameters [b epsilon alpha beta alpha_c beta_c]

pbounds = [0 0 0 0.05 0 0.05;
    1 1 1 13 1 13];    % bounds [lower; upper] * parameters [b epsilon alpha beta alpha_c beta_c]


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

% model names, and corresponding parameters and subplot label
names = {'RW model' 'RW model', 'RR model', 'WSLS model',...
    'RW-CK', 'RW-CK', 'RW-CK', 'RW-CK'};
symbols = {'\alpha' '\beta' 'b' '\epsilon' '\alpha' '\beta' '\alpha_c' '\beta_c'};
subplotLabel = {'A' 'B' 'C' 'D' 'E'};

%% Section 2a: simulate and fit data for RW

% simulate and fit data for nRep repetitions
for count = 1:nRep

    % random values for generating data
    alpha = rand;
    beta = exprnd(4);
    
    % simulate data
    [choice, reward, pt, q] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds, Npt);
    
    % fit the data
    [xf, LL] = fit_M3RescorlaWagner_v1(choice, reward, pt, pbounds(:, 3:4));

    % store true values, estimated values, and - loglikelihood
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
    [xf, LL] = fit_M1random_v1(choice, pbounds(:,1));
    
    % store true values, estimated values and the - loglikelihood
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
    [xf, LL] = fit_M2WSLS_v1(choice, reward, pbounds(:,2));
    
    % store true values, estimated values and the - loglikelihood
    fminX.sim(4,count) = epsilon;
    fminX.fit(4,count) = xf;
    fminX.negLL(3,count) = LL;
    
    clear xf LL choice reward
end

%% Section 2d: simulate and fit data for RW-CK model

for count = 1:nRep
    
    % random values for generating data
    alpha = rand;
    beta = exprnd(4);
    alpha_c = rand;
    beta_c  = exprnd(4);
    
    % simulate data
    [choice, reward, pt, q] = simulate_M4RWCK_v1(T, alpha, beta, alpha_c,...
        beta_c, rprob, rbounds, Npt);
    
    % fit the data
    [xf, LL] = fit_M4RWCK_v1(choice, reward, pt, pbounds(:, 3:6));

    % store true values, estimated values, and - loglikelihood
    fminX.sim(5,count) = alpha;
    fminX.sim(6,count) = beta;
    fminX.sim(7,count) = alpha_c;
    fminX.sim(8,count) = beta_c;
    fminX.fit(5,count) = xf(1);
    fminX.fit(6,count) = xf(2);
    fminX.fit(7,count) = xf(3);
    fminX.fit(8,count) = xf(4);
    fminX.negLL(4,count) = LL;

    % clear repeating variables from the workspace
    clear xf LL choice reward pt
end


%% Section 3: basic parameter recovery plots (true vs estimated)

% initiate figure
fh = figure('Name', 'Parameter Recovery'); clf;
set(gcf, 'Position', [300   313   1100   650])

% plot simulate versus fit parameters
for i = 1:size(fminX.sim,1)
    % select the subplot
    subplot(2,ceil(size(fminX.sim,1)/2),i); hold on
    
    % plot the true vs fit parameter values
    plot(fminX.sim(i,:), fminX.fit(i,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
    
    % plot the diagonal
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
end

% mark the 'bad' parameter values
for i = 1:size(fminX.sim,1)
    % determine threshold for 'bad' parameter values
    if i == 1 || i == 3 || i == 4
        thresh = 0.25;
    elseif i == 2 || i == 6 || i == 8
        thresh = 5;
    end
    
    % parameter values that exceed the threshold
    ind = abs(fminX.sim(i,:) - fminX.fit(i,:)) > thresh;
    
    % select the subplot
    subplot(2,ceil(size(fminX.sim,1)/2),i); hold on
    
    % set 'bad' parameter values in grey
    plot(fminX.sim(i,ind), fminX.fit(i,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
    'markerfacecolor', [1 1 1]*0.5)
end

% set titles and axis labels
for i = 1:size(names, 2)
    % select subplot
    subplot(2,ceil(size(fminX.sim,1)/2),i); hold on;
    
    % select the parameter name
    t(i) = title(names{i});
    
    % add label
    xlabel(sprintf('true %s', symbols{i}))
    ylabel(sprintf('estimated %s', symbols{i}));

    % set font and tick sizes
    ax = gca;
    set(ax, 'tickdir', 'out', 'fontsize', 12);
    
    % add subplot label
  %  addABCs(ax, [-0.06 0.05], 12, subplotLabel(i));
end

% save figure
if savePlots
    fh.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'Parameter_recovery', 'fit-vs-true.png');
    saveas(gcf, filename)
end

%% Section 4: Check correlation between the estimated alpha and beta parameters

% RW
% create table containing parameter values
alpha = fminX.fit(1,:)';
beta = fminX.fit(2,:)';
parTable = table(alpha, beta);

% plot the table and check if the correlation is significant
corrplot(parTable, 'type', 'Spearman', 'testR','on', 'rows', 'pairwise');
% th = findall(fh_cor(1).Parent.Parent, 'type', 'text', 'String', '{\bf Correlation Matrix}'); 
% th.String = 'RW';
%annotation('textbox', [.02 .01 .29 .05], 'String', 'estimated parmeters');

% save figure
if savePlots
%     fh2.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'Parameter_recovery', 'fit-vs-fit_RW.png');
    saveas(gcf, filename)
end
% ---

% RW CK
% create table containing parameter values
alpha = fminX.fit(5,:)';
beta = fminX.fit(6,:)';
alphaCK = fminX.fit(7,:)';
betaCK = fminX.fit(8,:)';
parTable = table(alpha, beta, alphaCK, betaCK);

% plot the table and check if the correlation is significant

corrplot(parTable, 'type', 'Spearman', 'testR','on', 'rows', 'pairwise');
% th = findall(fh_cor(1).Parent.Parent, 'type', 'text', 'String', '{\bf Correlation Matrix}'); 
% th.String = 'RW';
%annotation('textbox', [.02 .01 .29 .05], 'String', 'estimated parmeters');

% save figure
if savePlots
%     fh2.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'Parameter_recovery', 'fit-vs-fit_RW.png');
    saveas(gcf, filename)
end

%% Section 5: Plot bias-variance tradeoff in parameter estimation (RMSE)

% calculate the error between the true and estimated parameter values
for i = 1:length(symbols)
    
    % calculate the error between the estimated and true parameter
    % values
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

 % initiate figure
 fh3 = figure('Name', 'Parameter Estimation Error'); clf;
 set(gcf, 'Position',[500   113   1050   370])

 % plot estimation error
 for i = 1:size(fminX.error,1)
     % select the subplot
     subplot(1,size(fminX.error,1),i); hold on;
     
     % boxplot
     boxplot(fminX.error(i,:)', 'Labels',{[]}); hold on;
     
     % set y axis limits
     yl = get(gca, 'ylim');
     yl_max = max(abs(yl));
     ylim([-yl_max yl_max]);
    
     % add reference line
     line([0, 3], [0, 0], 'LineStyle', ':', 'Color', 'black');
     
     % add label
     xlabel(sprintf('%s', symbols{i}));
     ylabel(sprintf('true - estimated error'));
     
     % set font and tick sizes
     ax = gca;
     set(ax, 'tickdir', 'out', 'fontsize', 11);
     
     % add subplot label
%      addABCs(ax, [-0.06 0.05], 11, subplotLabel(i));
 end

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Parameter_recovery', 'biasVariance_allModels.png');
    saveas(gcf, filename)
end

% done