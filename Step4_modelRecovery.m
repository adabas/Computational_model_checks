% The script assesses (in terms of probability) the extent to which each
% model fits to the data simulated from every model. This allows us to
% arbitrate between the different models of interest. 
% This will be summarised in a confusion matrix. The confusion matix is
% defined as the probability that simulated by one model is best fit by
% another, i.e., p(fit model|simulated model).
%
% This is a modification of the script by Rob Williams and Anne Collins
% (2019).
%
% Modified by Aroma Dabas
%   - version       February 2020
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% set seed
rng(2, 'twister');

% ================== Modify ===============================================

T    = 132;          % number of trials
mu   = [0.7 0.4];    % reward probabilities
nRep = 100;          % number of repetitions#
rbounds = [0 1];     % bounds of the mean reward

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath('./HelperFunctions')
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots = 1;
plotFolder = './Figures/ModelSimulation/';

%% Section 2: Confusion matrix.
% simulate and fit data. Calculate the fit probability.

% initiate confusion matrix 
CM = zeros(3);

for count = 1:nRep
    count
    
    % open figure
    figure(1); clf;
    
    % calculate probability
    FM = round(100*CM/sum(CM(1,:)))/100;
    
    % display number/text on scaled image 
    t = imageTextMatrix(FM);
    
    % set font color of values less than 0.3 to white
    set(t(FM'<0.3), 'color', 'w')
    set(t, 'fontsize', 22)
    
    hold on;
    
    % add matrix lines
    [l1, l2] = addFacetLines(CM);
    
    % add count number as title
    title(['count = ' num2str(count)]);
   
    % set ticks and labels
    set(gca, 'xtick', [1:3], 'ytick', [1:3], 'fontsize', 28, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')
    
    % update the figure
    drawnow
    
    % Model 1
    b = rand;
    [a, r, pt] = simulate_M1random_v1(T, rbounds, b, mu);
    [BIC, iBEST, BEST] = fit_all_v1(a, r, pt);
    CM(1,:) = CM(1,:) + BEST;
    
    % Model 2
    epsilon = rand;
    [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, mu);
    [BIC, iBEST, BEST] = fit_all_v1(a, r, pt);
    CM(2,:) = CM(2,:) + BEST;
    
    % Model 3
    alpha = rand;
    beta  = 1 + exprnd(1);
    [a, r, pt] = simulate_M3RescorlaWagner_v1(T, alpha, beta, mu, rbounds);
    [BIC, iBEST, BEST] = fit_all_v1(a, r, pt);
    CM(3,:) = CM(3,:) + BEST;
    
end

% settings
title(sprintf('confusion matrix:\np(fit model | simulated model)'))
set(gcf, 'Position', [311   217   700   600])
set(gca, 'fontsize', 24);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Model_recovery', 'confusionMatrix.png');
    saveas(gcf, filename)
end

%% Section 3: Inversion matrix
% Given the model that fits our data bestm which model is likely to have
% generated the data? p(simulated data|fit model).

for i = 1:size(CM,2)
    iCM(:,i) = CM(:,i) / sum(CM(:,i));
end

% open figure
figure(2); clf
set(gcf, 'Position', [211   17   700   600]);

 % display number/text on scaled image 
t = imageTextMatrix(round(iCM, 2));

% set font color of values less than 0.3 to white
set(t(iCM<0.3), 'color', 'w')
set(t, 'fontsize', 22)

% add matrix lines
hold on
[l1, l2] = addFacetLines(CM);

 % set ticks and labels
set(gca, 'xtick', [1:3], 'ytick', [1:3], 'fontsize', 28, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('simulated model')
title(sprintf('inversion matrix:\np(simulated model | fit model)'))
set(gca, 'fontsize', 24);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Model_recovery', 'inversionMatrix.png');
    saveas(gcf, filename)
end
