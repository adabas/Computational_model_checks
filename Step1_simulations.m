% The script simulates data for the specified models, and plots
% model-independent measures (p(correct), trial-by-trial choices and p(stay))
%
% Explanation of the parameters for the 5 models:
%       Model 1: Random responding
%               b : random choice at mid point
%       Model 2: Win-stay-lose-shift
%               epsilon : chooses the option probabilistically that is
%               rewarded and switches away from unrewarded
%       Model 3: Rescorla Wagner
%               alpha : learning rate
%               beta : inverse temperature function
%       Model 5: Choice kernel (people tend to repeat their previous actions)
%               alpha_c : learning rate for the choice kernel; choice
%               kernel keeps track of how frequently an option was choosen
%               in the recent past
%               beta_c : inverse temperature function for choice kernel
%       Model 4: Rescorla-Wagner + choice kernel
%               alpha
%               beta
%               alpha_c
%               beta_c
%
% =========================================================================

clear; close all

%% Section 1: Settings (modify this section)

% set up colors
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

% add paths
addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath(genpath('./HelperFunctions'))

% set seed for result reproducibility
rng(7, 'twister');

% experiment parameters
T       = 96;       % number of trials
rbounds = [0 1];    % bounds of the mean reward
nrep    = 50;       % number of simulation repetitions
rprob   = [0.8 0.8 0.3 0.3]; % reward probability for [HR LR] stimuli 
Npt     = 0;        % number of partial trials

% specify the number and name of the model that you are interested in
% simulating
nModels     = [1 2 3 4 5];
nameModels  = {'Random responding', 'Noisy WSLS', 'RW', 'RW-CK', 'CK'};

% for plotting
savePlots   = 0;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

%% Section 2: simulate

% Model 1: Random responding
for n = 1:nrep
    b = 0.5;    % initiate bias parameter
    [a, r] = simulate_M1random(T, rbounds, b, rprob);
    sim(1).a(:,n) = a;  
    sim(1).r(:,n) = r;
end
clear a r

% Model 2: Noisy Win-stay-lose-shift
for n = 1:nrep
    epsilon = 0.3; % probability of selecting an option
    [a, r] = simulate_M2WSLS(T, rbounds, epsilon, rprob); 
    sim(2).a(:,n) = a;
    sim(2).r(:,n) = r;
end
clear a r

% Model 3: Rescorla Wagner
for n = 1:nrep
    alpha   = 0.1; 
    beta    = 5; 
    [a, r] = simulate_M3RescorlaWagner(T, alpha, beta, rprob, rbounds);
    sim(3).a(:,n) = a;
    sim(3).r(:,n) = r;
end
clear a r

% Model 4: Rescorla Wagner with Choice kernel
for n = 1:nrep
    alpha   = 0.1; 
    beta    = 5;
    alpha_c = 0.1;
    beta_c = 3;
    [a, r] = simulate_M4RWCK(T, alpha, beta, alpha_c, beta_c, rprob, rbounds);
    sim(4).a(:,n) = a;
    sim(4).r(:,n) = r;
end
clear a r

% Model 5: Choice kernel
for n = 1:nrep
    alpha_c = 0.1;
    beta_c = 3;
    [a, r] = simulate_M5CK(T, alpha_c, beta_c, rprob, rbounds);
    sim(5).a(:,n) = a;
    sim(5).r(:,n) = r;
end
clear a r

%% Section 2: Plot choice selection for HR stimuli
% Is the choice behaviour evolving over the trials to move towards the HR
% threshold? Or is the choice behaviour completely random for such a task?

% ======================= Plot choice probability =========================

% convert choices (coded 1 to 4) to reward conditions: HR: 1 & LR: 0
for i_mod = 1:length(nModels)
    sim(i_mod).a(sim(i_mod).a == 2) = 1;
    sim(i_mod).a(sim(i_mod).a == 3 | sim(i_mod).a == 4) = 0;
end

% initiate the plot
fh.score = figure('Name','score'); 
set(fh.score,'position',[10 90 1400 400], 'paperunits','centimeters','Color','w');
set(gca, 'fontsize', 12)

% loop over the number of models 
for m = 1:length(nModels)
    
    hold on
    
    % open the model's subplot
    subplot(1, length(nModels), m)
    
    % store choices for the model
    highRewChoice = sim(m).a; %choice = sim(m).a;
        
    % calculate the probability that HR was selected
    score = (nansum(highRewChoice == 1, 1))/T;
    
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
        title(sprintf('Model %i: %s \n alpha = %.2f;\n beta = %.2f', m, nameModels{m}, alpha, beta));
    elseif m == 4
        title(sprintf('Model %i: %s \n alpha = %.2f; beta = %.2f \n alpha_c = %.2f; beta_c = %.2f',...
            m, nameModels{m}, alpha, beta, alpha_c, beta_c));
    elseif m == 5
        title(sprintf('Model %i: %s \n alpha_c = %.2f; \n beta_c = %.2f',...
            m, nameModels{m}, alpha_c, beta_c));
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

% cycle over models
for m = 1:length(nModels)
    
    hold on
    % open the model's subplot
    subplot(ceil(length(nModels)/2), 2, m)
    
    % smoothing over the raw data
    smoothingkernel = 6;
    
    % calculate HR choice mean over each trial
    HR_choiceMean = nanmean(sim(m).a, 2);
    
    % add data to the plot
    line([0, length(HR_choiceMean)], [rprob(1), rprob(1)],...
        'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
    line([0, length(HR_choiceMean)], [rprob(2), rprob(2)],...
        'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
    plot(HR_choiceMean, ':', 'color', AZred, 'linewidth',0.5)
    plot(mySmooth(HR_choiceMean', smoothingkernel,[], 'backward'),...
        '-','color', AZred,'linewidth',1); 

    % y axis limits
    ylim([-0.1 1.1]);
    
    % add labels
    ylabel('p(HR stimulus)');
    xlabel('trial');
    legend({'p(reward|HR)', 'p(reward|LR)', 'mean choice', 'smoothed mean choice'}, 'location', 'southeast')
    
    if m == 1
        title(sprintf('Model %i: %s \n b = %.2f', m, nameModels{m}, b));
    elseif m == 2
        title(sprintf('Model %i: %s \n epsilon = %.2f', m, nameModels{m}, epsilon));
    elseif m == 3
        title(sprintf('Model %i: %s \n alpha = %.2f; beta = %.2f', m, nameModels{m}, alpha, beta));
    elseif m == 4
        title(sprintf('Model %i: %s \n alpha = %.2f; beta = %.2f \n alpha_c = %.2f; beta_c = %.2f',...
            m, nameModels{m}, alpha, beta, alpha_c, beta_c));
    elseif m == 5
        title(sprintf('Model %i: %s \n alpha_c = %.2f; beta_c = %.2f',...
            m, nameModels{m}, alpha_c, beta_c));
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

%% p(stay)
% determine the times the agent reselects the option that previously gave a
% pleasant, & unpleasant reward. Here we use the choice conditions instead
% of the individual 4 choices.

% loop over the simulations
for i = 1:length(sim)
    
    % loop over repetitions
    for n = 1:nrep
        
        sim(i).wsls(:,n) = analysis_WSLS(sim(i).a(:,n)', sim(i).r(:,n)');
        
    end
    
    % take mean of the wsls probability over all repetitions
    wsls(:,i) = nanmean(sim(i).wsls,2);
end

%% Plot WSLS (condition version) behavior for all models

% open figure
figure(3); clf; hold on;

% plot
l = plot(wsls);

% adjust limits
ylim([0 1])

% set plot and tick properties
set(l, 'marker', '.', 'markersize', 50, 'linewidth', 3)
set(gca, 'xtick', [1 2], 'XTickLabel', {'Unpleasant', 'Pleasant'}, ...
    'tickdir', 'out', 'fontsize', 14, 'xlim', [0.5 3.5])

% set title and labels
legend(nameModels, 'Location', 'southeast')
xlabel('previous event type')
ylabel('p(stay)')

% save plot
if savePlots
    fh.data.PaperPositionMode = 'auto';
    filename = fullfile(plotFolder, sprintf('./WSLS/WSLS_models_%i-%i.png', nModels(1), nModels(end)));
    saveas(gcf, filename)
end

% DONE