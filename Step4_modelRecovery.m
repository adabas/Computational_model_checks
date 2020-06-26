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

T    = 100;         % number of trials
mu   = [0.7 0.4];   % reward probabilities
nRep = 100;         % number of repetitions
rbounds = [0 1];    % bounds of the mean reward
Npt  = 0;           % number of partial trials
nMod = 3;           % number of models
% type = 1;           % 1 for determining best fitting model using BIC. 0 for
%                     % using negative log likelihood.
pbounds = [0 0 0 0;
    1 1 1 7];    % bounds [lower; upper] * parameters [b epsilon alpha beta]
   

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath(genpath('./HelperFunctions'))
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots = 1;
plotFolder = './Figures/ModelSimulation/';

% colors
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;
plotCol = {AZred AZcactus AZsky};

%% Section 2a: Confusion matrix.
% simulate and fit data. Calculate the best-fit probability conditional on the model used for simulation.

% initiate confusion matrix 
CM = zeros(nMod);

% initiate the plot
fh.cm = figure('Name','CM');
set(fh.cm,'position',[100 50 700 600],'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

for count = 1:nRep
    
    % Model 1
    b = rand;
    [a, r, pt] = simulate_M1random_v1(T, rbounds, b, mu, Npt);
    [~, ~, BEST, pars] = fit_all_v1(a, r, pt, nMod, pbounds);
    CM(1,:) = CM(1,:) + BEST;
    
    % Model 2
    epsilon = rand;
    [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, mu, Npt);
    [~, ~, BEST] = fit_all_v1(a, r, pt, nMod, pbounds);
    CM(2,:) = CM(2,:) + BEST;
    
    % Model 3
    alpha = rand;
    beta  = 5 + exprnd(5); % 1 + exprnd(1);
    [a, r, pt] = simulate_M3RescorlaWagner_v1(T, alpha, beta, mu, rbounds, Npt);
    [~, ~, BEST] = fit_all_v1(a, r, pt, nMod, pbounds);
    CM(3,:) = CM(3,:) + BEST;
    
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
    set(gca, 'xtick', [1:3], 'ytick', [1:3], 'fontsize', 18, ...
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
% estimation are too similar or if they are actually very dissimilar and
% the results you plotted thus far are more likely related to the 
% parameter estimation procedure. 

% It might alse be useful to plot the learning performance similarly to
% Step1_simulations for the current parameter settings and models. If all
% models show a similar performance, then it is likely that the parameter
% range is not good, because we know that these models should in principle
% make clearly dissociable predictions. 

% Thanks for the suggestions. I tried the suggestions in section 2b and 3.

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
    filename = fullfile(plotFolder, 'Model_recovery', 'CM.png');
    saveas(gcf, filename)
end

%% Section 2b: Additional analysis
% Simulate choices, and estimate data using constricted parameters.
% Use the log likelihood values for the best fitting model to create the CM.
% I determined the parameter bounds for eliciting similar learning performance
% by trial and testing the Step1_simulation.m script.

% -- Poor learning behaviour --
% specify the parameter bounds
pbounds = [0.4 0.6 0.6 0.5;
    0.6 1 1 2];

% name the figure
nameFig = 'CM (low learning)';

% plot the confusion matrix
CM_plot(T, mu, 'name', nameFig, 'pbounds', pbounds, 'nMod', nMod,...
    'nRep', nRep);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Model_recovery', 'CM_lowLearning.png');
    saveas(gcf, filename)
end

% -- Fast learning behaviour --
% specify the parameter bounds.
pbounds = [0.3 0.1 0.7 4;
    0.4 0.4 0.8 6];

% name the figure
nameFig = 'CM (fast learning)';

% plot the confusion matrix
CM_plot(T, mu, 'name', nameFig, 'pbounds', pbounds, 'nMod', nMod,...
    'nRep', nRep);

% save plot
if savePlots
    filename = fullfile(plotFolder, 'Model_recovery', 'CM_fastLearning.png');
    saveas(gcf, filename)
end


%% Section 3: Plot learning performance
% Simulate data using a model and estimate each model's learning
% performance. Are the model orignial and estimation-based results disparate from one another? 

% parameter bounds for simulating the data
pbounds = [0.3 0.1 0.02 5;  % 0.7 4
    0.4 0.4 0.2 30]; %  0.8 20

% initiate the plot
fh.ch = figure('Name','Model based choice estimation');
set(fh.ch,'position',[10 100 800 1000],'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

% plot
for model = 1:nMod
    
    % simulate and estimate data
    [sim, fit] = choiceEstimation('pbounds', pbounds, 'model', model, 'rewProb', mu,...
        'nMod', nMod, 'ntrials', T);
    
    % plot
    % smoothing over the raw data
    smoothingkernel = 6;

    % extract high-reward choice option
    highRewAction = find(mu==max(mu));

    hold on
    % open the model's subplot
    subplot(nMod-1, 2, model)

    % calculate mean HR choice over each trial for simulated data
    HR_choiceMean = nanmean(sim(1).a, 2);

    % rescale the choice mean such that it ranges between 0 (LR) and 1 (HR)
    if highRewAction == 2
        HR_choiceMean = (2 - HR_choiceMean)';
    else
        HR_choiceMean = (HR_choiceMean - 1)';
    end

    % calculate mean HR choice over each trial for fitted data
    for m = 1:nMod
        % calculate mean
       fit(m).mean = nanmean(fit(m).a,2);

       % rescale
       if highRewAction == 2
            fit(m).mean = (2 - fit(m).mean)';
       else
            fit(m).mean = (fit(m).mean - 1)';
       end

    end

    % add data to the plot
    plot(mySmooth(HR_choiceMean, smoothingkernel,[], 'backward'),...
        '-','color', plotCol{model},'linewidth',2); hold on
    for m = 1:nMod
        plot(mySmooth(fit(m).mean, smoothingkernel, [], 'backward'),...
            '-', 'color', plotCol{m}, 'linewidth', 0.8); hold on
    end

     % y axis limits
    ylim([-0.1 1.1]);

    % add labels
    ylabel('p(HR stimulus)');
    xlabel('trial');
    legend({sprintf('sim. model %i', model), 'Model 1', 'Model 2', 'Model 3'}, 'location', 'southeast')
    
    legend boxoff

end

if savePlots
    filename = fullfile(plotFolder, 'Model_recovery', 'choice_estimated.png');
    saveas(gcf, filename)
end

% For the given parameter bounds, the models' learning performance is quite
% similar, particularly for data simulated by model 2.

%% Section 4: Inversion matrix
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
    filename = fullfile(plotFolder, 'Model_recovery', 'IM.png');
    saveas(gcf, filename)
end

% DONE