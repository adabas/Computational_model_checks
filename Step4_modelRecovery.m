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
% pbounds = [0 0.05 0 0.05 0;       % parameter bounds updated to empirical data     
%   1 1 200 1 40];
pbounds = [0 0.05 0 0.05 0;       % parameter bounds updated to empirical data     
  1 1 50 1 50];


% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
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

%% Section 2a: Confusion matrix.
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
    [a, r, s] = simulate_M2WSLS_v2(T, rbounds, epsilon, mu);
    [~, BESTbic, ~, BESTnegll, ~] = fit_all_v2(a, r, nMod, pbounds, s);
   % [~, ~, BEST] = fit_all_v2(a, r, nMod, pbounds, s);
    CMbic(1,:) = CMbic(1,:) + BESTbic;
    CMnll(1,:) = CMnll(1,:) + BESTnegll;
    
    % Model 3
    alpha = rand;
    beta  = 3 + exprnd(3);
    [a, r, s] = simulate_M3RescorlaWagner_v2(T, alpha, beta, mu, rbounds);
    [~, BESTbic, ~, BESTnegll, ~] = fit_all_v2(a, r, nMod, pbounds, s);
   % [~, ~, BEST] = fit_all_v2(a, r, pt, nMod, pbounds, s);
    CMbic(2,:) = CMbic(2,:) + BESTbic;
    CMnll(2,:) = CMnll(2,:) + BESTnegll;
    
     % Model 4
    alpha = rand;
    beta  = 3 + exprnd(3); 
    alpha_c = rand;
    beta_c  = 3 + exprnd(3); 
    [a, r, s] = simulate_M4RWCK_v2(T, alpha, beta, alpha_c, beta_c, mu, rbounds);
    [~, BESTbic, ~, BESTnegll, ~]  = fit_all_v2(a, r, nMod, pbounds, s);
    %[~, ~, BEST] = fit_all_v2(a, r, pt, nMod, pbounds, s); 
    CMbic(3,:) = CMbic(3,:) + BESTbic;
    CMnll(3,:) = CMnll(3,:) + BESTnegll;
    
    % Model 5
    alpha_c = rand;
    beta_c  = 3 + exprnd(3);
    [a, r, s] = simulate_M5CK_v2(T, alpha_c, beta_c, mu, rbounds);
    [~, BESTbic, ~, BESTnegll, ~] =  fit_all_v2(a, r, nMod, pbounds, s);
    %[~, ~, BEST] = fit_all_v2(a, r, pt, nMod, pbounds, s);
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

%% Section 2b: Additional analysis
% Simulate choices, and estimate data using constricted parameters.
% Use the log likelihood values for the best fitting model to create the CM.
% I determined the parameter bounds for eliciting similar learning performance
% by trial and testing the Step1_simulation.m script.

% -- Poor learning behaviour --
% % specify the parameter bounds
% pbounds = [0 0 0.7 4 0.7 4;
%      1 1 0.8 6 0.8 6];
% 
% % name the figure
% nameFig = 'CM (low learning)';
% 
% % plot the confusion matrix
% CM_plot(T, mu, 'name', nameFig, 'pbounds', pbounds, 'nMod', nMod,...
%     'nRep', nRep);
% 
% % save plot
% if savePlots
%     filename = fullfile(plotFolder, 'Model_recovery', 'CM_lowLearning.png');
%     saveas(gcf, filename)
% end
% 
% % -- Fast learning behaviour --
% % specify the parameter bounds.
% pbounds = [0.3 0.1 0.7 4 0.7 4;
%     0.4 0.4 0.8 6 0.8 6];
% 
% % name the figure
% nameFig = 'CM (fast learning)';
% 
% % plot the confusion matrix
% CM_plot(T, mu, 'name', nameFig, 'pbounds', pbounds, 'nMod', nMod,...
%     'nRep', nRep);
% 
% % save plot
% if savePlots
%     filename = fullfile(plotFolder, 'Model_recovery', 'CM_fastLearning.png');
%     saveas(gcf, filename)
% end


%% Section 3: Plot learning performance (but we already plotted this in Step1)
% Simulate data using a model and estimate each model's learning
% performance. Are the model original and estimation-based results different from one another? 

% % parameter bounds for simulating the data
% pbounds = [0 0 0.7 4 0.7 4;
%      1 1 0.8 6 0.8 6];

% % initiate the plot
% fh.ch = figure('Name','Model based choice estimation');
% set(fh.ch,'position',[10 100 800 1000],'paperunits','centimeters','Color','w');
% set(gca, 'fontsize', 12)
% 
% % plot
% for model = 1:nMod
%     
%     % simulate and estimate data
%     [sim, fit] = choiceEstimation('pbounds', pbounds, 'model', model, 'rewProb', mu,...
%         'nMod', nMod, 'ntrials', T);
%     
%     % plot
%     % smoothing over the raw data
%     smoothingkernel = 6;
% 
%     % extract high-reward choice option
%     highRewAction = find(mu==max(mu));
% 
%     hold on
%     % open the model's subplot
%     subplot(nMod-1, 2, model)
% 
%     % calculate mean HR choice over each trial for simulated data
%     HR_choiceMean = nanmean(sim(1).a, 2);
% 
%     % rescale the choice mean such that it ranges between 0 (LR) and 1 (HR)
%     if highRewAction == 2
%         HR_choiceMean = (2 - HR_choiceMean)';
%     else
%         HR_choiceMean = (HR_choiceMean - 1)';
%     end
% 
%     % calculate mean HR choice over each trial for fitted data
%     for m = 1:nMod
%         % calculate mean
%        fit(m).mean = nanmean(fit(m).a,2);
% 
%        % rescale
%        if highRewAction == 2
%             fit(m).mean = (2 - fit(m).mean)';
%        else
%             fit(m).mean = (fit(m).mean - 1)';
%        end
% 
%     end
% 
%     % add data to the plot
%     plot(mySmooth(HR_choiceMean, smoothingkernel,[], 'backward'),...
%         '-','color', plotCol{model},'linewidth',2); hold on
%     for m = 1:nMod
%         plot(mySmooth(fit(m).mean, smoothingkernel, [], 'backward'),...
%             '-', 'color', plotCol{m}, 'linewidth', 0.8); hold on
%     end
% 
%      % y axis limits
%     ylim([0.3 0.9]);
% 
%     % add labels
%     ylabel('p(HR stimulus)');
%     xlabel('trial');
%     legend({sprintf('sim. model %i', model), 'Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5'}, 'location', 'southeast')
%     
%     legend boxoff
% 
% end
% 
% if savePlots
%     filename = fullfile(plotFolder, 'Model_recovery', 'choice_estimated.png');
%     saveas(gcf, filename)
% end

%% Section 4: Inversion matrix
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