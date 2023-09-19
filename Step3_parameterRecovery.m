% The script tests parameter recovery for the models: simple RW, WSLS, RW-CK and CK. This plots:
% 1. Scatter plot with fitted vs simulated parameters
% 2. Correlation between parameters of a) RW model and b) RW-CK model
% 3. bias-variance trade off: is there is an overall bias in the estimated
%    parameters, and if the variance of the estimated parameters is widely
%    spread.
% 
% This code is a modification of the code provided by Bob Wilson and Anne
% Collins [2019].
%
% Updates by Aroma Dabas [dabas@cbs.mpg.de]
%   01-2023     Cleaned and removed random responding model
%   01-2020     Version 1
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% set seed
rng(234, 'twister');

% ================== Modify ===============================================
% experiment parameters
T       = 96;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nRep    = 50;       % number of simulation repetitions
rprob   = [0.8 0.8 0.3 0.3]; % reward probability for [HR LR] stimuli
Npt     = 0;        % number of partial trials
pbounds = [0.5 0 0.05 0 0.05 0;     % parameter bounds updated to empirical data     
  0.5 1 1 20 1 20];          % fixing null models parameter space to 0.5

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

savePlots   = 0;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

% model names, and corresponding parameters and subplot label
names = {'RW model' 'RW model', 'WSLS model','RW-CK', 'RW-CK', 'RW-CK',...
    'RW-CK', 'CK', 'CK'};
symbols = {'\alpha' '\beta' '\epsilon' '\alpha' '\beta' '\alpha_c' '\beta_c',...
    '\alpha_c', '\beta_c'};
subplotLabel = {'A' 'B' 'C' 'D' 'E'};

%% Section 2a: simulate and fit data for RW

% simulate and fit data for nRep repetitions
for count = 1:nRep

    % random values for generating data
    alpha = rand;
    beta = exprnd(4);
    
    % simulate data
    [choice, reward, s] = simulate_M3RescorlaWagner_v2(T, alpha, beta, rprob, rbounds);
    
    % fit the data
    [xf, LL] = fit_M3RescorlaWagner_v2(choice, reward, pbounds(:, 3:4), s);

    % store true values, estimated values, and - loglikelihood
    fminX.sim(1,count) = alpha;
    fminX.sim(2,count) = beta;
    fminX.fit(1,count) = xf(1);
    fminX.fit(2,count) = xf(2);
    fminX.negLL(count) = LL;

    % clear repeating variables from the workspace
    clear xf LL choice reward

end

%% Section 2b: simulate and fit data for noisy WSLS model

for count = 1:nRep
    
    % random parameter value for generating data
    epsilon = rand;
    
    % simulate data
    [choice, reward, s] = simulate_M2WSLS_v2(T, rbounds, epsilon, rprob);
    
    % fit the data
    [xf, LL] = fit_M2WSLS_v2(choice, reward, pbounds(:,2), s);
    
    % store true values, estimated values and the - loglikelihood
    fminX.sim(3,count) = epsilon;
    fminX.fit(3,count) = xf;
    fminX.negLL(3,count) = LL;
    
    clear xf LL choice reward
end

%% Section 2c: simulate and fit data for RW-CK model

for count = 1:nRep
    
    % random values for generating data
    alpha = rand;
    beta = exprnd(4);
    alpha_c = rand;
    beta_c  = exprnd(4);
    
    % simulate data
    [choice, reward, s] = simulate_M4RWCK_v2(T, alpha, beta, alpha_c, beta_c,...
        rprob, rbounds);
    
    % fit the data
    [xf, LL] = fit_M4RWCK_v2(choice, reward, pbounds(:, 3:6), s);

    % store true values, estimated values, and - loglikelihood
    fminX.sim(4,count) = alpha;
    fminX.sim(5,count) = beta;
    fminX.sim(6,count) = alpha_c;
    fminX.sim(7,count) = beta_c;
    fminX.fit(4,count) = xf(1);
    fminX.fit(5,count) = xf(2);
    fminX.fit(6,count) = xf(3);
    fminX.fit(7,count) = xf(4);
    fminX.negLL(4,count) = LL;

    % clear repeating variables from the workspace
    clear xf LL choice reward pt
end

%% Section 2d: simulate and fit data for CK model

for count = 1:nRep
    
    % random values for generating data
    alpha_c = rand;
    beta_c  = exprnd(4);
    
    % simulate data
    [choice, reward, s] = simulate_M5CK_v2(T, alpha_c, beta_c, rprob, rbounds);
    
    % fit the data
    [xf, LL] = fit_M5ChoiceKernel_v2(choice, reward, pbounds(:, 5:6), s);

    % store true values, estimated values, and - loglikelihood
    fminX.sim(8,count) = alpha_c;
    fminX.sim(9,count) = beta_c;
    fminX.fit(8,count) = xf(1);
    fminX.fit(9,count) = xf(2);
    fminX.negLL(5,count) = LL;

    % clear repeating variables from the workspace
    clear xf LL choice reward pt
end

%% Section 3: parameter recovery correlation (true vs estimated)

% estimate Kendall's tau
for iparam = 1:numel(names)
    [rho(iparam), pval(iparam)] = corr(fminX.sim(iparam,:)', fminX.fit(iparam,:)',...
        'Type', 'Kendall');
end
fprintf('\nKendall tau correlation between simulated and fit parameters:\n');
fprintf('    RW model:   alpha rho = %.2f (p = %.2f), beta rho = %.2f (p = %.2f)\n', rho(1), pval(1), rho(2), pval(2));
fprintf('    WSLS model: epsilon rho = %.2f (p = %.2f)\n', rho(3), pval(3));
fprintf('    CK model:   alpha rho = %.2f (p = %.2f), beta rho = %.2f (p = %.2f)\n', rho(8), pval(8), rho(9), pval(9));
fprintf('    RW-CK model:\n');
fprintf('        RW:     alpha rho = %.2f (p = %.2f), beta rho = %.2f (p = %.2f)\n', rho(4), pval(4), rho(5), pval(5));
fprintf('        CK:     alpha rho = %.2f (p = %.2f), beta rho = %.2f (p = %.2f)\n', rho(6), pval(6), rho(7), pval(7));

% --- PLOT
% initiate figure
fh = figure('Name', 'Parameter Recovery'); clf;
set(gcf, 'Position', [300   313   1100   650])

% plot simulate versus fit parameters
for i = 1:size(fminX.sim,1)
    % select the subplot
    subplot(2,ceil(size(fminX.sim,1)/2),i); hold on
    
    % plot the true vs fit parameter values
    plot(fminX.sim(i,:), fminX.fit(i,:), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1)
    
    % fit linear regression line
    h1 = lsline;
    h1.Color = 'r';   
%     % plot the diagonal
%     xl = get(gca, 'xlim');
%     plot(xl, xl, 'k--')
end

% mark the 'bad' parameter values
for i = 1:size(fminX.sim,1)
    % determine threshold for 'bad' parameter values
    if i == 1 || i == 3 || i == 4 || i == 6 || i == 8
        thresh = 0.25;
    elseif i == 2 || i == 5 || i == 7 || i == 9
        thresh = 10;
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
end

% save figure
if savePlots
    fh.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'ParameterRecovery', 'fit-vs-true.png');
    saveas(gcf, filename)
end

%% Section 4: Check correlation between the estimated alpha and beta parameters

% RW
% create table containing parameter values
alpha = fminX.fit(1,:)';
beta = fminX.fit(2,:)';
parTable = table(alpha, beta);

% plot the table and check if the correlation is significant
figure;
corrplot(parTable, 'type', 'Spearman', 'testR','on', 'rows', 'pairwise');
title('RW model correlation plot')

% save figure
if savePlots
%     fh2.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'ParameterRecovery', 'RW_fittedAlphaBetaCorrelation.png');
    saveas(gcf, filename)
end
clear partable
% ---

% RW CK
% create table containing parameter values
alpha = fminX.fit(4,:)';
beta = fminX.fit(5,:)';
alphaCK = fminX.fit(6,:)';
betaCK = fminX.fit(7,:)';
parTable = table(alpha, beta, alphaCK, betaCK);

% plot the table and check if the correlation is significant
figure;
corrplot(parTable, 'type', 'Spearman', 'testR','on', 'rows', 'pairwise');
title('RW-CK model correlation plot')
% save figure
if savePlots
%     fh2.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, 'ParameterRecovery', 'RWCK_fittedAlphaBetaCorrelation.png');
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
fprintf('\nRMSE values of RW model:\nalpha = %.3f\n', fminX.RMSE(1,:))
fprintf('beta = %.3f\n', fminX.RMSE(2,:))

% ---
% plot 1: alpha of all models

% convert all alpha values to long format
error.alpha = [fminX.error(1,:)'; fminX.error(4,:)';...
    fminX.error(6,:)'; fminX.error(8,:)'];

name.alpha = cell(200, 1);  % Initialize the cell array
for i = 1:nRep; name.alpha{i} = 'RW'; end
for i = nRep+1:nRep*2; name.alpha{i} = 'RW-CK (RW)'; end
for i = nRep*2+1:nRep*3; name.alpha{i} = 'RW-CK (CK)'; end
for i = nRep*3+1:nRep*4; name.alpha{i} = 'CK'; end

% initiate figure
fh_alpha = figure('Name', 'Alpha parameter recovery', 'color', 'w');
ax = axes('NextPlot','add','FontSize',16,'TickDir','out');

% boxplot
boxplot(error.alpha, name.alpha, 'Parent', ax) %, 'Symbol', 'o', 'Whisker', inf)

h =findobj(ax, 'LineStyle','--'); set(h, 'LineStyle','-');
h = findobj('Marker','+'); set(h, 'Marker','o'); set(h,'MarkerFaceColor','b', 'MarkerEdgeColor', 'b')

% modify boxes
myboxes = findobj(ax,'Tag','Box');
arrayfun( @(box) patch( box.XData, box.YData, 'b', 'FaceAlpha', 0.5), myboxes(1:4) )

if savePlots
    filename = fullfile(plotFolder, 'ParameterRecovery', 'biasVariance_alpha.png');
    saveas(gcf, filename)
end

% ---
% plot 2: beta of all models
% convert all alpha values to long format
error.beta = [fminX.error(2,:)'; fminX.error(5,:)';...
    fminX.error(7,:)'; fminX.error(9,:)'];

name.beta = cell(200, 1);  % Initialize the cell array
for i = 1:nRep; name.beta{i} = 'RW'; end
for i = nRep+1:nRep*2; name.beta{i} = 'RW-CK (RW)'; end
for i = nRep*2+1:nRep*3; name.beta{i} = 'RW-CK (CK)'; end
for i = nRep*3+1:nRep*4; name.beta{i} = 'CK'; end

% initiate figure
fh_beta = figure('Name', 'Beta parameter recovery', 'color', 'w');
ax = axes('NextPlot','add','FontSize',16,'TickDir','out');

% boxplot
boxplot(error.beta, name.beta, 'Parent', ax) %, 'Symbol', 'o', 'Whisker', inf)

h =findobj(ax, 'LineStyle','--'); set(h, 'LineStyle','-');
h = findobj('Marker','+'); set(h, 'Marker','o'); set(h,'MarkerFaceColor','b', 'MarkerEdgeColor', 'b')

% modify boxes
myboxes = findobj(ax,'Tag','Box');
arrayfun( @(box) patch( box.XData, box.YData, 'b', 'FaceAlpha', 0.5), myboxes(1:4) )

if savePlots
    filename = fullfile(plotFolder, 'ParameterRecovery', 'biasVariance_beta.png');
    saveas(gcf, filename)
end

% plot 3: epsilon WSLS model

error.epsilon = fminX.error(3,:)';

name.epsilong = cell(50, 1);  % Initialize the cell array
for i = 1:nRep; name.epsilon{i} = 'WSLS'; end

% initiate figure
fh_epsilon = figure('Name', 'Epsilon parameter recovery', 'color', 'w');
ax = axes('NextPlot','add','FontSize',16,'TickDir','out');

% boxplot
boxplot(error.epsilon, name.epsilon, 'Parent', ax) %, 'Symbol', 'o', 'Whisker', inf)

h =findobj(ax, 'LineStyle','--'); set(h, 'LineStyle','-');
h = findobj('Marker','+'); set(h, 'Marker','o'); set(h,'MarkerFaceColor','b', 'MarkerEdgeColor', 'b')

% modify boxes
myboxes = findobj(ax,'Tag','Box');
arrayfun( @(box) patch( box.XData, box.YData, 'b', 'FaceAlpha', 0.5), myboxes(1) )
ylim([-.4 .4])

if savePlots
    filename = fullfile(plotFolder, 'ParameterRecovery', 'biasVariance_epsilon.png');
    saveas(gcf, filename)
end

% done