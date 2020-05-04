% The script assesses (in terms of probability) the extent to which each
% model fits to the data simulated from every model. This allows us to
% arbitrate between the different models of interest. 
% This will be summarised in a confusion matrix. The confusion matix is
% defined as the probability that simulated by one model is best fit by
% another, i.e., p(fit model|simulated model).
%
% This is a modification of the script by Bob Wilson and Anne Collins
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
Npt  = 32;          % number of partial trials

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath(genpath('./HelperFunctions'))
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots = 1;
plotFolder = './Figures/ModelSimulation/';

%% Section 2: Confusion matrix.
% simulate and fit data. Calculate the best-fit probability conditional on the model used for simulation.

% initiate confusion matrix 
CM = zeros(3);

for count = 1:nRep
    
    % Model 1
    b = rand;
    [a, r, pt] = simulate_M1random_v1(T, rbounds, b, mu, Npt);
    [~, ~, BEST] = fit_all_v1(a, r, pt);
    CM(1,:) = CM(1,:) + BEST;
    
    % Model 2
    epsilon = rand;
    [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, mu, Npt);
    [~, ~, BEST] = fit_all_v1(a, r, pt);
    CM(2,:) = CM(2,:) + BEST;
    
    % Model 3
    alpha = rand;
    beta  = 1 + exprnd(1);
    [a, r, pt] = simulate_M3RescorlaWagner_v1(T, alpha, beta, mu, rbounds, Npt);
    [~, ~, BEST] = fit_all_v1(a, r, pt);
    CM(3,:) = CM(3,:) + BEST;
    
    % I changed the order of simulation/estimation and plotting. If I
    % understand it correctly, you previously didn't plot the last
    % iteration (100) but instead iteration 0. Now you plot from 1 - 100.
    
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
    
end

% This is a nice way of checking model recovery. I would certainly keep on
% using this method but would also suggest an extension. 
% It looks like the probability of the best-fitting model is relatively
% similar when you simulate the data with model 3. Give that this is the
% most important/interesting model, I would argue that we should 
% dig into this a little bit deeper.

% One potential explanation is that the model likelihoods for the three
% models are very similar. Another explanation is that the likelihoods are 
% actually sufficiently different but the parameter estimation is not
% sufficiently accurate, which is not super likely, because you already 
% conducted the parameter recovery analysis. A third explanation is that 
% the parameter of the simulations are in a range where they cannot be
% reliably dissociated. 

% To figure out which explanation is most likely, I would suggest the
% following: Conduct an additional analysis where you don't estimate the
% free parameters but where you use the same parameters as for the
% simulation to compute the log likelihood. Then you use the same confusion
% matrix but not with the best BIC but the best log-likelihood. That way we have a
% better idea if the model likelihoods independent of the parameter
% estimation are too similar or if they are actually very similar and
% the results you plotted thus far are more likely related to the 
% parameter estimation procedure. 

% It might alse be useful to plot the learning performance similarly to
% Step1_simulations for the current parameter settings and models. If all
% models show a similar performance, then it is likely that the parameter
% range is not good, because we know that these models should in principle
% make clearly dissociable predictions. 

% In the long run, I would also include exceedance probabilites, i.e.,
% simply checking which model is most likely on the group level. I'd to
% this both for model recovery with simulations and later for the empirical data.
% But let's first concentrate on the first extension I proposed. 

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
% Given the model that fits our data best, which model is most likely to have
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
