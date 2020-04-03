% Fitting LSim data to the Rescorla Wagner model and Softmax function.
% Ensure you are in the parent directory of the script.
% 
% The learning model used to fit data.choices is a simple Rescorla-Wagner
% (Rescorla & Wagner 1972) function, using a softmax
% data.choice/link/observation function:  
%   learning model (RW):     vA <- vA + alpha*(r-vA)
%   function (softmax):    p(A) = exp(beta*vA)/(exp(beta*vA)_exp(beta*vB))
%
% Code by Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 12/03/19
% =========================================================================
% 
% Still to consider/do:
%   1. what to do with the data.choices that the participant did not make but
%   the system selected because the participants were not fast enough?
% 
% Version 1:    15-03-19    Plots only for individual data
% Version 2:                Plot for group analysis
% Version 3:    25-03-19    Moved group analysis code to
%                           LSim_plotData_group.m
% Version 4:    30-10-19    Included alpha and beta likelihood surface maps
% Version 5:    23-01-20    Normalised pleasantness ratings
% Version 6:    28-01-20    Included fmincon optimiser for parameter
%                           estimation.
%                           Updated grid plot to include fmincon parameter
%                           estimates.
%                           Supperimposed estimated trial-by-trial
%                           probability of selecting HR stimuli using the
%                           fmincon best parameter estimates.

%% Section 1: Preparation

close all;
clearvars

rng(244, 'twister');    % set seed

% ================== Modify ===============================================
subjects    = 0186;     % specify subject ID
savePlots   = false;    % true will save plots in plotFolder directory
plotFolder  = "./Figures/SubjectLevel";     % figure path as a string

% specify the bounds, and the number of bins for grid search
bounds  = [0 1; % alpha
    0 50]; % beta
n.bin   = [20 30] ;

% store labels for plotting
cond        = {'High Reward', 'Low Reward'};    % unconditioned stimuli labels
n.cond      = length(cond);                     % no. of UCS
paramLabel  = {'alpha', 'beta'};                % parameter labels
n.param     = length(paramLabel);               % no. of parameters 

% ================== Add paths ============================================

% add path of the current folder
tmp = fileparts(which('Step5_plotData_SubjectLevel'));
addpath(tmp);

% add path to required folders in the current folder
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'LikelihoodFunctions'))
addpath(fullfile(tmp, 'SimulationFunctions'))

% add path to the data folder
rootdir     = tmp(1:end-32);
datapath    = fullfile(rootdir, '01_Data');

% ================== Load subject's data ==================================

% add path to the subjects's data folder
addpath(fullfile(datapath, sprintf('subject_0%i', subjects)));

% load 'results' matrix
load(sprintf('subject_0%i_sessionII_results.mat',subjects), 'results');

% load settings
load(sprintf('subject_0%i_settingsS2.mat',subjects));

% load stimuli
load(sprintf('subject_0%i_sessionII_items.mat', subjects), 'itemmat');

% ================== Plot colors ==========================================
global AZred AZcactus AZsky

AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

%% Section 2: Data restructuring and organisation

% store data.choices: HR coded as 1 and LR coded as 2
data.choice = cell2mat(results(2:end,9));

% store the number of trials
n.trials = settings.design.ntrialblock * settings.design.nblock;

% calculate the probability with which the participant selected HR stimuli
data.score = (nansum(data.choice == 1))/(n.trials);

% load the pleasantness ratings
for i = 1:n.trials
    if isempty(results{i+1,14})
       ratings(i) = NaN;    % replace the partial trial events with NaN
    else
       ratings(i) = results{i+1,14};
    end
end
data.rate.raw = ratings';

% normalise the ratings
data.rate.min    = min(data.rate.raw);     % store the minimum    
data.rate.max    = max(data.rate.raw);     % store the maximum
data.rate.norm   = (data.rate.raw - data.rate.min)/...
    (data.rate.max - data.rate.min);

% store the event types shown on each trial (Positive, Negative-Neutral & 'Partial')
data.eventType   = results(2:end, 12);
idPartial  = find(cellfun('isempty', data.eventType));  % the empty cells are partial trials
data.eventType(idPartial) = {'Partial'};

% store the liking ratings of the six selected stimuli
for i = 1:length(itemmat)
   v(i) = cell2mat(itemmat(i).liking);
end

% store the corresponding condition ids of the six stimuli: HR coded as 1,
% LR coded as 2 and baseline coded as 3
for i = 1:length(itemmat)
    cond_id(i,1) = itemmat(i).cond_id;
end

% average the liking for the HR and LR conditions
data.v0(1) = mean(v(cond_id == 1));  % HR
data.v0(2) = mean(v(cond_id == 2));  % LR

% clear temporary variables
clear cond_id v i

%% Section 3: Plot choice behaviour
% Is the data.choice behaviour evolving over the trials to move towards the 
% p(reward|HR) threshold?

% ================== Plot p(correct data.choices)  =========================

% plot using the barScatter.m function
fh.score = figure('Name','score');set(fh.score,'position',[100 500 250 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');
barScatter(data.score,[],[],true);   % barplot function with error bars
set(gca,'xtick',[]);
ylabel('p(HR)');
title(sprintf('Subject %i', subjects));
ylim([0 1])
box off  

% save figure
if savePlots
    fh.score.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/prob_0%i.png', plotFolder, subjects))
end

% ======== Plot raw behavioural data & smoothed response function ========

data.choice_plot = 2-data.choice(:)';     % recode data.choices such that HR is stored 
                                % as 1 and LR is stored as 0
smoothingkernel = 6;            % no. of trials over which to smooth

% start plot
fh.choice = figure('Name', 'Trial-by-trial choice');
box off; hold on;
set(fh.choice,'position', [500 500 700 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');

% set y axis limits
ylim([-0.1 1.1]);

% add data to the plot
line([0, length(data.choice_plot)],...
    [settings.design.rerate{1}, settings.design.rerate{1}],...
    'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
line([0, length(data.choice_plot)],...
    [settings.design.rerate{2},settings.design.rerate{2}],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
plot(data.choice_plot, 'k*', 'markersize',5);
plot(mySmooth(data.choice_plot,smoothingkernel,[],'backward'),...
        '-','color', [0.8500, 0.3250, 0.0980],'linewidth',1); 

% add labels
title(sprintf('subject %i',subjects));
ylabel('p(HR stimuli)');
xlabel('trial');
legend({'p(reward|HR)', 'p(reward|LR)','choice','smoothed choice (HR)'},'location','northeastoutside')
legend boxoff

if savePlots
    fh.choice.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/choices_0%i.png', plotFolder, subjects))
end

clear data.choice_plot smoothingkernel

%% Section 4:
% Find the parameter values of RW+softmax function that maximise the
% likelihood of the data. As a sanity check, run both fmincon and grid
% search methods

if size(bounds,2)~=n.param
    error('number of bounds and parameters don''t match')
end

% Run fmincon optimiser
% =====================================================================
% set fmincon options: 
    % - 'MaxFunEval': maximum number of function evaluations allowed.
    % - 'Display' set to 'notify': displays ouptput if the function does
    %                               not converge.
    % - 'algorithm' set to 'active-set': can take large steps, which adds
    %                               speed. The algorithm is effective on
    %                               some problems with nonsmooth constraints.
options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% run optimization over 10 starting points
for count = 1:10

    % random starting points [alpha beta]
    X0  = [rand exprnd(10)];

    % create RW+softmax function handle to input in fmincon
    obFunc = @(x) lik_M3RescorlaWagner_v1(data.choice', data.rate.norm', x(1), x(2), idPartial);
%         obFunc = @(x) lik_M3RescorlaWagner_v1(data.choice', data.rate.norm, X0(1), X0(2), idPartial);

    % store the lower and upper bounds of [alpha beta] parameters
    LB = [bounds(1,1) bounds(2,1)];
    UB = [bounds(1,2) bounds(2,2)];

    % run fmincon to check the best fitting parameter to the data
    [xf, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

    % store starting value (Xint), fitted value (Xfit) and the negative
    % LL values (negLL)
    fminX.initiate(1,count) = X0(1);
    fminX.initiate(2,count) = X0(2);
    fminX.fit(1,count)      = xf(1);
    fminX.fit(2,count)      = xf(2);
    fminX.negLL(count)      = NegLL;

    % clear repeating variables from the workspace
    clear X0 xf NegLL LB UB

end

% find global best
[mf,i]=min(fminX.negLL(:));
fminX.pars = fminX.fit(:,i);


% Run grid search method
% =====================================================================

% create bins for both the parameters
for iParam = 1:n.param
    range = linspace(bounds(iParam,1),bounds(iParam,2),n.bin(iParam)+1);
    p{iParam} = range(2:end); % stay just off the zero bounds
end

% select each bin and determine the likelihood of the data given the
% parameters in that bin
params = nan(1,n.param);
for t = 1:(n.bin(1))
    params(1) = p{1}(t);    % alpha value
    for tt = 1:(n.bin(2))
        params(2) = p{2}(tt);   % beta value

        % determine the log likelihood value for the parameter
        llh(t,tt) = -lik_M3RescorlaWagner_v1(data.choice',...
            data.rate.norm', params(1), params(2), idPartial);

        %[llh2(t,tt)] = RL_fitmodel(params, data.choice,...
        %    data.rate.norm, idPartial, [0.5], [0.5]);

    end
end


%% Section 5: Plot estimated parameters


% Plot grid surface with estimated parameters from both grid search and fmincon method 
% =====================================================================

fh.grid = figure('Name','Grid');
set(fh.grid,'position',[600 50 650 650],'paperunits','centimeters','Color','w');

% find minimum and maximum values from the grid search
mi      = min(llh(:));
[ma,i]  = max(llh(:));

% create matrices containing bins of alpha and beta parameters
x=repmat(1:length(p{1,1}),length(p{1,2}),1)';   % repeating 
y=repmat(1:length(p{1,2}),length(p{1,1}),1);

% plot grid surface
imagesc(p{1,1}(1:end),p{1,2}(1:end),llh',[mi,ma])
colorbar

% plot the best fitting parameters from the grid search
hold on
plot(p{1,1}(x(i)), p{1,2}(y(i)), 'ok')

% plot the best fitting parameters from the fmincon optimiser function
plot(fminX.pars(1),fminX.pars(2),'*k')

% add labels and edit settings
xlabel('alpha')
ylabel('beta')
set(gca,'fontsize',14)
title(sprintf('Subject %i',subjects))

dim = [0.65 0 0.5 0.05];
annotation('textbox',dim, 'String', '*: fmincon; o: grid search',...
    'EdgeColor', 'none', 'FontSize', 11);

% save plot
if savePlots
    fh.grid.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/GridPlots/M3_RescorlaWagner+Softmax/0%i.png', plotFolder, subjects))
end

% remove unnecessary variables
clear x y t tt i

% Plot estimated choice probability using the best fitting parameters
% =====================================================================

% estimate the trial-by-trial probability using the best parameter
% estimates obtained from the fmincon method.

[~, PP, d, Q] = lik_M3RescorlaWagner_v1(data.choice',...
            data.rate.norm', fminX.pars(1), fminX.pars(2), idPartial);

% select the trial-by-trial choices plot
figure(fh.choice)

% plot the estimated trial-by-trial choices
plot(PP(:,1)','--','color', [0, 0.4470, 0.7410],'linewidth',0.7);

% update legend
legend({'p(reward|HR)', 'p(reward|LR)','choice','smoothed choice (HR)', 'estimated choice (HR)'},...
    'location','northeastoutside')
legend boxoff

% save plot
if savePlots
    fh.choice.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/M3_RescorlaWagner+Softmax/Est_choices_0%i.png', plotFolder, subjects))
end

% Plot prediction error
% =========================================================================

fh.delta = figure('Name', 'Trial-by-trial prediction error');
box off; hold on;
set(fh.delta,'position', [300 300 700 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');

% set y axis limits
ylim([-1 1]);

% add data to the plot
line([0, length(data.choice_plot)],...
    [0, 0],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.3);  hold on    
plot(d,...
        '-','color', AZsky,'linewidth',1);

% add labels
title(sprintf('subject %i',subjects));
ylabel('prediction error');
xlabel('trial');

if savePlots
    fh.delta.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/M3_RescorlaWagner+Softmax/PredictionError_0%i.png', plotFolder, subjects))
end


% DONE