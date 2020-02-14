% The script simulates data for the specified models, and plots
% model-independent measures (p(correct), trial-by-trial choices and p(stay))
%
% Explanation of the parameters for the 5 models:
%       Model 1: Random responding
%               b : bias for an option; here it is assumed that
%               participants choose the options randomly, perhaps with some
%               overall bias for one option (PE or NE stimuli) over the other.
%               Reward value ranges from 0 (negative reward) to 1 (positive
%               reward) with 0.5 as neutral reward.
%       Model 2: Win-stay-lose-shift
%               epsilon : chooses the option probabilistically that is
%               rewarded and switches away from unrewarded
%       Model 3: Rescorla Wagner
%               alpha : learning rate
%               beta : inverse temperature function
%       Not implemented yet, but upcoming:
%       Model 4: Choice kernel (people tend to repeat their previous actions)
%               alpha_c : learning rate for the choice kernel; choice
%               kernel keeps track of how frequently an option was choosen
%               in the recent past
%               beta_c : inverse temperature function for choice kernel
%       Model 5: Rescorla-Wagner + choice kernel
%               alpha
%               beta
%               alpha_c
%               beta_c
%       Model 6: Rescorla Wagner + stickiness function
%
% Updated by Aroma Dabas [dabas@cbs.mpg.de]
%   12-2019         First version
%   01-2020         Implemented random responding model
%                   Implemented Noisy WSLS model
%                   Commented and cleaned the scripts
% 
% Original code by:
% Bob Wilson & Anne Collins
% 2018
% Code to produce figure 2 in submitted paper "Ten simple rules for the
% computational modeling of behavioral data"

clear; close all

%% Section 1: Settings (modify this section)

% set up colors
global AZred AZcactus AZsky

AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

% add paths
addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath('./HelperFunctions')

% set seed for result reproducibility
rng(7, 'twister');

% experiment parameters
T       = 132;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nrep    = 40;       % number of simulation repetitions
rprob   = [0.7 0.4]; % reward probability for [HR LR] stimuli 

% specify the number and name of the model that you are interested in
% simulating.
nModels     = [1 2 3];
nameModels  = {'Random responding', 'Noisy WSLS', 'RW+Softmax'};

% for plotting
savePlots   = 1;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

%% Section 2: simulate
% Model 1: Random responding
for n = 1:nrep
    b = 0.5;    % initiate bias parameter
    [a, r] = simulate_M1random_v1(T, rbounds, b, rprob);    % simulate data
    sim(1).a(:,n) = a;  
    sim(1).r(:,n) = r;
end
clear a r

% Model 2: Noisy Win-stay-lose-shift
for n = 1:nrep
    epsilon = 0.1; % probability of selecting an option
    [a, r] = simulate_M2WSLS_v1(T, rbounds, epsilon, rprob);
    sim(2).a(:,n) = a;
    sim(2).r(:,n) = r;
end
clear a r

% Model 3: Rescorla Wagner
for n = 1:nrep
    alpha   = 0.2;
    beta    = 12;
    [a, r, ~, PP] = simulate_M3RescorlaWagner_v1(T, alpha, beta, rprob, rbounds);
    sim(3).a(:,n) = a;
    sim(3).r(:,n) = r;
end
clear a r

% % Model 4: Choice kernel
% for n = 1:Nrep
%     alpha_c = 0.1;
%     beta_c = 3;
%     [a, r] = simulate_M4ChoiceKernel_v1(T, mu, alpha_c, beta_c);
%     sim(4).a(:,n) = a;
%     sim(4).r(:,n) = r;
% end
% clear a r
% 
% % Model 5: Rescorla-Wagner + choice kernel
% for n = 1:Nrep
%     alpha = 0.05;
%     beta = 4;
%     alpha_c = 0.1;
%     beta_c = 1;
%     [a, r] = simulate_M5RWCK_v1(T, mu, alpha, beta, alpha_c, beta_c);
%     sim(5).a(:,n) = a;
%     sim(5).r(:,n) = r;
% end
% clear a r

%% Section 2: Plot choice selection for HR stimuli
% Is the choice behaviour evolving over the trials to move towards the HR
% threshold? Or is the choice behaviour completely random for such a task?

% ======================= Plot choice probability =========================
% initiate the plot
fh.score = figure('Name','score'); 
set(fh.score,'position',[10 90 800 400], 'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

% loop over the number of models 
for m = 1:length(nModels)
    
    hold on
    % open the model's subplot
    subplot(1, length(nModels), m)
    
    % store choices for the model
    choice = sim(m).a;
    
    % convert the scoring such as the HR is 1 and LR is 0
    choice = 2 - choice; 

    % calculate the probability that HR was selected
    score = (nansum(choice == 1, 1))/T;
    
    % calculate overall probability that the high reward was selected
    barScatter(score,[],[],true);
    ylim([0 1])
    
    % some settings
    set(gca,'xtick',[]);
    ylabel('p(HR)', 'FontWeight', 'normal');
    if m == 1
        title(sprintf('Model %i: %s \n b = %.2f', m, nameModels{m}, b));
    elseif m == 2
        title(sprintf('Model %i: %s \n epsilon = %.2f', m, nameModels{m}, epsilon));
    elseif m == 3
        title(sprintf('Model %i: %s \n alpha = %.2f; beta = %.2f', m, nameModels{m}, alpha, beta));
    end

    box off  
    
end

% add footer
AxesH = axes('Units', 'Normalized', 'Position', [0,0,1,1], 'Visible', 'off', ...
             'NextPlot', 'add');
text(0.2, 0.07, sprintf('Simulated for n = %i', n), 'Parent', AxesH, ...
     'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

 % save plot
if savePlots
    fh.score.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, sprintf('./Choices/prob_models_%i-%i.png', nModels(1), nModels(end)));
    saveas(gcf, filename)
end

% ================== Plot smoothed trial wise choice ======================

% initiate the plot
fh.data = figure('Name','Data');
set(fh.data,'position',[10 100 800 800],'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

for m = 1:length(nModels)
    
    hold on
    % open the model's subplot
    subplot(length(nModels)-1, 2, m)
    
    % smoothing over the raw data
    smoothingkernel = 6;
    
    % calculate choice mean over each trial
    choiceMean = nanmean(sim(m).a, 2);
    
    % rescale the choice mean such that it ranges between 0 (LR) and 1 (HR)
    choiceMean = (2 - choiceMean)';
    
    % add data to the plot
    line([0, length(choiceMean)], [rprob(1), rprob(1)],...
        'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
    line([0, length(choiceMean)], [rprob(2), rprob(2)],...
        'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
    plot(choiceMean, ':', 'color', AZred,'linewidth',0.5)
    plot(mySmooth(choiceMean,smoothingkernel,[],'backward'),...
        '-','color', AZred,'linewidth',1); 

    % y axis limits
    ylim([-0.1 1.1]);
    
    % add labels
    ylabel('p(HR stimuli)');
    xlabel('trial');
    legend({'p(reward|HR)', 'p(reward|LR)', 'mean choice', 'smoothed mean choice'},'location','southeast')
    
    if m == 1
        title(sprintf('Model %i: %s \n b = %.2f', m, nameModels{m}, b));
    elseif m == 2
        title(sprintf('Model %i: %s \n epsilon = %.2f', m, nameModels{m}, epsilon));
    elseif m == 3
        title(sprintf('Model %i: %s \n alpha = %.2f; beta = %.2f', m, nameModels{m}, alpha, beta));
    end
    
    legend boxoff
    
end

% add footer
AxesH = axes('Units', 'Normalized', 'Position', [0,0,1,1], 'Visible', 'off', ...
             'NextPlot', 'add');
text(0.2, 0.05, sprintf('Simulated for n = %i', n), 'Parent', AxesH, ...
     'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% save plot
if savePlots
    fh.data.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, sprintf('./Choices/choices_models_%i-%i.png', nModels(1), nModels(end)));
    saveas(gcf, filename)
end


%% p(stay), also referred to as win stay lose shift (WSLS) analysis
% determine the times the agent reselects the option that previously gave a
% pleasant, neutral & unpleasant reward.

% loop over the simulations
for i = 1:length(sim)
    
    % loop over repetitions
    for n = 1:nrep
        
        sim(i).wsls(:,n) = analysis_WSLS_v1(sim(i).a(:,n)', sim(i).r(:,n)');
        
    end
    
    % take mean of the wsls probability over all repetitions
    wsls(:,i) = nanmean(sim(i).wsls,2);
end

%% Plot WSLS behavior for all models

% open figure
figure(3); clf; hold on;

% plot
l = plot(wsls);

% adjust limits
ylim([0 1])

% set plot and tick properties
set(l, 'marker', '.', 'markersize', 50, 'linewidth', 3)
set(gca, 'xtick', [1 2 3], 'XTickLabel', {'Unpleasant', 'Neutral', 'Pleasant'}, ...
    'tickdir', 'out', 'fontsize', 14, 'xlim', [0.5 3.5])

% set title and labels
legend({'M1: Random responding' 'M2: Noisy WSLS' 'M3: RW+softmax'}, ...
    'location', 'southeast')
xlabel('previous event type')
ylabel('p(stay)')

% save plot
if savePlots
    fh.data.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, sprintf('./WSLS/WSLS_models_%i-%i.png', nModels(1), nModels(end)));
    saveas(gcf, filename)
end

% DONE