% Fitting data from LSim3 study to the RW model.
% 
% Code by Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 12/03/19
% =========================================================================
%% Section 1: Preparation

close all;
%clearvars

rng(244, 'twister');    % set seed

% ================== Modify ===============================================
subjects    = 60;     % specify subject ID
savePlots   = true;    % true will save plots in plotFolder directory
saveData    = true;
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
addpath(fullfile(tmp, 'HelperFunctions'));
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'LikelihoodFunctions'))
addpath(fullfile(tmp, 'SimulationFunctions'))
addpath(genpath(fullfile(tmp, 'DataOutput')))

% add path to the data folder
rootdir     = tmp(1:end-32);
datapath    = fullfile(rootdir, '01_Data', 'LSim_3_behavioural');

% ================== Load subject's data ==================================
addpath(fullfile(datapath, sprintf('subject_00%i', subjects)));

% load 'results' matrix
load(sprintf('subject_%i_P3_results.mat',subjects), 'results');

% load settings
load(sprintf('subject_%i_settingsP1.mat',subjects));

% load stimuli
load(sprintf('subject_%i_sessionII_items.mat', subjects), 'itemmat');

% ================== Plot colors ==========================================
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

%% Section 2: Data restructuring and organisation

% store data.choices: HR coded as 1 and LR coded as 2
data.choice = cell2mat(results(2:end,9));

% store the number of trials
n.trials = settings.design.ntrialblock * settings.design.nblock;

% calculate the probability with which the participant selected HR stimuli
data.score = (sum(data.choice == 1))/(n.trials - sum(isnan(data.choice)));%(n.trials); %(nansum(data.choice == 1))/(n.trials);

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

% store the corresponding condition ids of the six stimuli: HR coded as 1,
% LR coded as 2 and baseline coded as 3
for i = 1:length(itemmat)
    cond_id(i,1) = itemmat(i).cond_id;
end

% clear temporary variables
clear cond_id i

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
    'paperposition',[0 0 5 5],'Color','w');

% set y axis limits
ylim([-0.1 1.1]);
%xlim([0, 100]);

% smoothed data
data.choice_smooth = mySmooth(data.choice_plot,smoothingkernel,[],'backward');

line([0, length(data.choice_plot)],...
    [settings.design.rerate{1}, settings.design.rerate{1}],...
    'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
line([0, length(data.choice_plot)],...
    [settings.design.rerate{2},settings.design.rerate{2}],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
plot(data.choice_plot, 'k*', 'markersize',5);
plot(data.choice_smooth,...
        '-','color', [0.8500, 0.3250, 0.0980],'linewidth',1);
line([(settings.design.ntrialblock+1), (settings.design.ntrialblock+1)],...
   [-0.1,1.1],'LineStyle', '-', 'Color', 'black', 'linewidth',0.5);  hold on
line([((settings.design.ntrialblock*2)+1), ((settings.design.ntrialblock*2)+1)],...
   [-0.1,1.1],'LineStyle', '-', 'Color', 'black', 'linewidth',0.5);  hold on
line([((settings.design.ntrialblock*3)+1), ((settings.design.ntrialblock*3)+1)],...
   [-0.1,1.1],'LineStyle', '-', 'Color', 'black', 'linewidth',0.5);  hold on

% add labels
title(sprintf('subject %i',subjects));
ylabel('p(HR stimuli)');
xlabel('trial');
legend({'p(reward|HR)', 'p(reward|LR)','choice','smoothed choice (HR)'},'location','southeast')
legend boxoff

if savePlots
    fh.choice.PaperPositionMode = 'auto';
    %print -depsc2 finalPlot1.eps
    saveas(gcf, sprintf('%s/Choices/choices_0%i.png', plotFolder, subjects))
end

clear data.choice_plot smoothingkernel

%% Save Data

% filenames
simCSV.fileName.probHR = fullfile('DataOutput', 'probabilityHR', 'simulationTask', 'simulation_pHR.csv');
simCSV.fileName.trialHR = fullfile('DataOutput', 'probabilityHR', 'simulationTask', 'simulation_trialwiseHR.csv');

% create data as tables
simCSV.table.probHR = table(subjects, data.score);
simCSV.table.probHR.Properties.VariableNames = {'Subjects', 'probHR'};
simCSV.table.trialHR = table(repelem(subjects, length(data.choice_plot))', (1:length(data.choice_plot))', data.choice_plot', data.choice_smooth');
simCSV.table.trialHR.Properties.VariableNames = {'Subjects', 'Trial', 'RawChoices', 'SmoothedChoices'};

if saveData
    % p(HR)
    if ~exist(simCSV.fileName.probHR, 'file')
        writetable(simCSV.table.probHR, simCSV.fileName.probHR);
    else
        simCSV.old.probHR = readtable(simCSV.fileName.probHR);
        if ~any(simCSV.old.probHR.Subjects == subjects)
            % concatenate table
            simCSV.new.probHR = [simCSV.old.probHR; simCSV.table.probHR];
            writetable(simCSV.new.probHR, simCSV.fileName.probHR);
        end
    end
    % trial-wise
    if ~exist(simCSV.fileName.trialHR, 'file')
        writetable(simCSV.table.trialHR, simCSV.fileName.trialHR);
    else
        simCSV.old.trialHR = readtable(simCSV.fileName.trialHR);
        if ~any(simCSV.old.trialHR.Subjects == subjects)
            % concatenate table
            simCSV.new.trialHR = [simCSV.old.trialHR; simCSV.table.trialHR];
            writetable(simCSV.new.trialHR, simCSV.fileName.trialHR);
        end
    end    
    
end

%% Section 4:
% Find the parameter values of RW+softmax function that maximise the
% likelihood of the data.

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

end

% find global best
[mf,i]=min(fminX.negLL(:));
fminX.pars = fminX.fit(:,i);

%% Section 5: Plot estimated parameters

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
    'location','southeast')
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
xlim([0, 100]);

% add data to the plot
line([0, length(data.choice_plot)], [0, 0],...
    'LineStyle', '--', 'Color', AZcactus, 'MarkerEdgeColor', 'none', 'linewidth',0.5);  hold on    
pe_line = plot(d,'-','color', AZsky,'MarkerEdgeColor', 'none', 'linewidth', 2.5);

% more modifications
set(gca, ...
    'FontName'   , 'Helvetica',...
    'TickDir'   , 'out' ,...
    'TickLength', [.01 .01],...
    'XColor'    , [.3 .3 .3],...
    'YColor'    , [.3 .3 .3],...
    'YTick'     , -0.8:0.2:0.8,...
    'XTick'     , 1:11:100,...
    'YGrid'     , 'on');

% add labels
%title(sprintf('subject %i',subjects));
ylabel('prediction error', 'FontSize', 15);
xlabel('trial', 'FontSize', 15);

if savePlots
    fh.delta.PaperPositionMode = 'auto';
    %print -depsc2 predictionError.eps
    %fh.delta.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/M3_RescorlaWagner+Softmax/PredictionError_0%i.png', plotFolder, subjects))
end

%% Section 6: Bandit task datasets

% load bandit task dataset
load(sprintf('subject_%i_bandit_results.mat',subjects), 'results');

% store data.choices: HR coded as 1 and LR coded as 2
data_band.choice = zeros((length(results)-1),1);
for i = 1:(length(results)-1)
    if strcmp(results{i+1,6}, 'NAN')
        data_band.choice(i) = NaN;
    else
        data_band.choice(i) = results{i+1,6};
    end
end
% store the number of trials
n.trials = length(data_band.choice);%settings.design.ntrialblock * settings.design.nblock;

% calculate the probability with which the participant selected HR stimuli
data_band.score = (nansum(data_band.choice == 1))/(n.trials);

% --- plot the choice behaviour ----

% ================== Plot p(correct data.choices)  =========================

% plot using the barScatter.m function
fh_bandit.score = figure('Name','bandit score');set(fh_bandit.score,'position',[100 500 250 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');
barScatter(data_band.score,[],[],true);   % barplot function with error bars
set(gca,'xtick',[]);
ylabel('p(HR)');
title(sprintf('Subject %i', subjects));
ylim([0 1])
box off  

% save figure
if savePlots
    fh_bandit.score.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/BanditTask/probBandit_0%i.png', plotFolder, subjects))
end

% ======== Plot raw behavioural data & smoothed response function ========

data_band.choice_plot = 2-data_band.choice(:)';     % recode data.choices such that HR is stored 
                                % as 1 and LR is stored as 0
smoothingkernel = 6;            % no. of trials over which to smooth

% start plot
fh_bandit.choice = figure('Name', 'Trial-by-trial choice (Bandit Task)');
box off; hold on;
set(fh_bandit.choice,'position', [500 500 700 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');

% set y axis limits
ylim([-0.1 1.1]);

%
data_band.choice_smooth = mySmooth(data_band.choice_plot,smoothingkernel,[],'backward');

% add data to the plot
line([0, length(data_band.choice_plot)],...
    [settings.design.rerate{1}, settings.design.rerate{1}],...
    'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
line([0, length(data_band.choice_plot)],...
    [settings.design.rerate{2},settings.design.rerate{2}],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
plot(data_band.choice_plot, 'k*', 'markersize',5);
plot(mySmooth(data_band.choice_plot,smoothingkernel,[],'backward'),...
        '-','color', [0.8500, 0.3250, 0.0980],'linewidth',1); 

% add labels
title(sprintf('subject %i',subjects));
ylabel('p(HR stimuli)');
xlabel('trial');
legend({'p(reward|HR)', 'p(reward|LR)','choice','smoothed choice (HR)'},'location','southeast')
legend boxoff

if savePlots
    fh_bandit.choice.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/BanditTask/choicesBandit_0%i.png', plotFolder, subjects))
end

%% dataset

% filenames
banditCSV.fileName.probHR = fullfile('DataOutput', 'probabilityHR', 'banditTask', 'bandit_pHR.csv');
banditCSV.fileName.trialHR = fullfile('DataOutput', 'probabilityHR', 'banditTask', 'bandit_trialwiseHR.csv');

% create data as tables
banditCSV.table.probHR = table(subjects, data_band.score);
banditCSV.table.probHR.Properties.VariableNames = {'Subjects', 'probHR'};
banditCSV.table.trialHR = table(repelem(subjects, length(data_band.choice_plot))', (1:length(data_band.choice_plot))', data_band.choice_plot', data_band.choice_smooth');
banditCSV.table.trialHR.Properties.VariableNames = {'Subjects', 'Trial', 'RawChoices', 'SmoothedChoices'};

if saveData
    % p(HR)
    if ~exist(banditCSV.fileName.probHR, 'file')
        writetable(banditCSV.table.probHR, banditCSV.fileName.probHR);
    else
        banditCSV.old.probHR = readtable(banditCSV.fileName.probHR);
        if ~any(banditCSV.old.probHR.Subjects == subjects)
            % concatenate table
            banditCSV.new.probHR = [banditCSV.old.probHR; banditCSV.table.probHR];
            writetable(banditCSV.new.probHR, banditCSV.fileName.probHR);
        end
    end
    
    % trial-wise
     if ~exist(banditCSV.fileName.trialHR, 'file')
        writetable(banditCSV.table.trialHR, banditCSV.fileName.trialHR);
    else
        banditCSV.old.trialHR = readtable(banditCSV.fileName.trialHR);
        if ~any(banditCSV.old.trialHR.Subjects == subjects)
            % concatenate table
            banditCSV.new.trialHR = [banditCSV.old.trialHR; banditCSV.table.trialHR];
            writetable(banditCSV.new.trialHR, banditCSV.fileName.trialHR);
        end
    end
end

% DONE