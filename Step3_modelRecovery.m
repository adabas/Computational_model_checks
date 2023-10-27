% The script tests model recovery. We test if we can arbitrate different
% models. This is summarised in a confusion matrix, p(fit model|simulated
% model). We also run the inverse of this matrix, i.e. inversion matrix
% p(simulated model | fit model).
%
% This is a modification of the script by Bob Wilson and Anne Collins
% (2019).
%
% Modified by Aroma Dabas
%   - version       February 2020
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% ================== Modify ===============================================

T       = 96;               % number of trials
mu      = [.8 .8 .3 .3];    % reward probabilities
nRep    =  50;              % number of repetitions
rbounds = [0 1];            % bounds of the mean reward
Npt     = 0;                % number of partial trials
nMod    = 4;                % number of models
pbounds = [0.5 0 0.05 0 0.05 0;     % parameter bounds updated to empirical data     
  0.5 1 1 20 1 20];    

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath(genpath('./HelperFunctions'))
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots = 0;
plotFolder = './Figures/ModelSimulation/';

% colors
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;
brown = [171 104 87]./255;
purple = [144 103 167]./255;
plotCol = {AZred AZcactus AZsky brown purple};

%% Section 2: Confusion matrix.
% simulate and fit data. Calculate the best-fit probability conditional on the model used for simulation.

% initiate confusion matrix using BIC
CMbic = zeros(nMod);
fh.cmbic = figure('Name','CM-BIC');
set(fh.cmbic,'position',[100 50 700 600],'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

% initiate confusion matrix using negative loglikelihood
CMnll = zeros(nMod);
fh.cmnll = figure('Name','CM-negLL');
set(fh.cmnll,'position',[100 50 700 600],'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

% set seed
rng(2, 'twister');

for count = 1:nRep
    
    % Model 2
    epsilon = rand;
    [a, r, s] = simulate_M2WSLS(T, rbounds, epsilon, mu);
    [~, ~, BESTbic, ~, BESTnegll] = fit_M2to5(a, r, nMod, pbounds, s);
    CMbic(1,:) = CMbic(1,:) + BESTbic;
    CMnll(1,:) = CMnll(1,:) + BESTnegll;
    
    % Model 3
    alpha = rand;
    beta  = 3 + exprnd(3);
    [a, r, s] = simulate_M3RescorlaWagner(T, alpha, beta, mu, rbounds);
    [~, ~, BESTbic, ~, BESTnegll] = fit_M2to5(a, r, nMod, pbounds, s);
    CMbic(2,:) = CMbic(2,:) + BESTbic;
    CMnll(2,:) = CMnll(2,:) + BESTnegll;
    
     % Model 4
    alpha = rand;
    beta  = 3 + exprnd(3); 
    alpha_c = rand;
    beta_c  = 3 + exprnd(3); 
    [a, r, s] = simulate_M4RWCK(T, alpha, beta, alpha_c, beta_c, mu, rbounds);
    [~, ~, BESTbic, ~, BESTnegll] = fit_M2to5(a, r, nMod, pbounds, s);
    CMbic(3,:) = CMbic(3,:) + BESTbic;
    CMnll(3,:) = CMnll(3,:) + BESTnegll;
    
    % Model 5
    alpha_c = rand;
    beta_c  = 3 + exprnd(3);
    [a, r, s] = simulate_M5CK(T, alpha_c, beta_c, mu, rbounds);
    [~, ~, BESTbic, ~, BESTnegll] = fit_M2to5(a, r, nMod, pbounds, s);
    CMbic(4,:) = CMbic(4,:) + BESTbic;
    CMnll(4,:) = CMnll(4,:) + BESTnegll;
    
    % ---
    % calculate probability for BIC
    figure(1);
    FMbic = round(100*CMbic/sum(CMbic(2,:)))/100;
    
    % display number/text on scaled image 
    t = imageTextMatrix(FMbic, {'WSLS', 'RW', 'RW-CK', 'CK'}, {'WSLS', 'RW', 'RW-CK', 'CK'});
    
    % set font color of values less than 0.3 to white
    set(t(FMbic'<0.3), 'color', 'w')
    set(t, 'fontsize', 22)
    
    hold on;
    
    % add matrix lines
    [l1, l2] = addFacetLines(CMbic);
    
    % add count number as title
    title(['count = ' num2str(count)]);
   
    % set ticks and labels
    set(gca, 'xtick', 1:nMod, 'ytick', 1:nMod, 'fontsize', 18, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')
    
    % update the figure
    drawnow
    
    % ---
    % calculate probability for negLL
    figure(2);
    FMnegLL = round(100*CMnll/sum(CMnll(2,:)))/100;
    
    % display number/text on scaled image 
    t = imageTextMatrix(FMnegLL, {'WSLS', 'RW', 'RW-CK', 'CK'}, {'WSLS', 'RW', 'RW-CK', 'CK'});
    
    % set font color of values less than 0.3 to white
    set(t(FMnegLL'<0.3), 'color', 'w')
    set(t, 'fontsize', 22)
    
    hold on;
    
    % add matrix lines
    [l1, l2] = addFacetLines(CMnll);
    
    % add count number as title
    title(['count = ' num2str(count)]);
   
    % set ticks and labels
    set(gca, 'xtick', 1:nMod, 'ytick', 1:nMod, 'fontsize', 18, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')
    
    % update the figure
    drawnow
end

% settings
figure(1)
title(sprintf('confusion matrix (BIC):\np(fit model | simulated model)'))
set(gcf, 'Position', [311   217   700   600]) 
set(gca, 'fontsize', 20);

figure(2)
title(sprintf('confusion matrix (-LL):\np(fit model | simulated model)'))
set(gcf, 'Position', [311   217   700   600]) 
set(gca, 'fontsize', 20);

% save plot
if savePlots
    figure(1)
    filename = fullfile(plotFolder, 'ModelRecovery', 'CM_BIC.png');
    saveas(gcf, filename);
    figure(2)
    filename = fullfile(plotFolder, 'ModelRecovery', 'CM_negLL.png');
    saveas(gcf, filename);
    
end

%% Section 3: Inversion matrix
% Given the model that fits our data best, which model is most likely to have
% generated the data? p(simulated data|fit model).

% ---
% BIC based IM
for i = 1:size(CMbic,2)
    iCMbic(:,i) = CMbic(:,i) / sum(CMbic(:,i));
end

% open figure
figure(3); clf
set(gcf, 'Position', [211   17   700   600]);

 % display number/text on scaled image 
t = imageTextMatrix(round(iCMbic, 2));

% set font color of values less than 0.3 to white
set(t(iCMbic<0.3), 'color', 'w')
set(t, 'fontsize', 20)

% add matrix lines
hold on
[l1, l2] = addFacetLines(iCMbic);

 % set ticks and labels
set(gca, 'xtick', [1:4], 'ytick', [1:4], 'fontsize', 28, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('simulated model')
title(sprintf('inversion matrix (BIC):\np(simulated model | fit model)'))
set(gca, 'fontsize', 22);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'ModelRecovery', 'IM_BIC.png');
    saveas(gcf, filename)
end

% ---
% negative LL based IM

for i = 1:size(CMnll,2)
    iCMnll(:,i) = CMnll(:,i) / sum(CMnll(:,i));
end

% open figure
figure(4); clf
set(gcf, 'Position', [211   17   700   600]);

 % display number/text on scaled image 
t = imageTextMatrix(round(iCMnll, 2));

% set font color of values less than 0.3 to white
set(t(iCMnll<0.3), 'color', 'w')
set(t, 'fontsize', 20)

% add matrix lines
hold on
[l1, l2] = addFacetLines(iCMnll);

 % set ticks and labels
set(gca, 'xtick', [1:4], 'ytick', [1:4], 'fontsize', 28, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('simulated model')
title(sprintf('inversion matrix (-LL):\np(simulated model | fit model)'))
set(gca, 'fontsize', 22);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'ModelRecovery', 'IM_negLL.png');
    saveas(gcf, filename)
end

% DONE