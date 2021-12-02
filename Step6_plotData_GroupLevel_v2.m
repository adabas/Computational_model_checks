% Fitting data to the following models.
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
%       Model 4: Rescorla-Wagner + choice kernel
%               alpha
%               beta
%               alpha_c
%               beta_c
%       Model 5: Choice Kernel
%               alpha_c
%               beta_c
% 
% Ensure you are in the parent directory of the script.
%
% Code by Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 04/2020
% =========================================================================

% NOTE: Check subjectAnalysis -> fit_all_v1
%% Section 1: Preparation
close all;
clearvars

rng(244);    % set seed

% ================== Modify ===============================================
subjects    = [11:18 20:21 23:55 57:60]; %[11:18, 20:21, 23:60];  % specify subject IDs
savePlots   = true;    % true will save plots in plotFolder directory
highRewAction = 2;      % set 2 and 1 to plot HR choices and LR choices, respectively.
rprob       = [0.8 0.3];
plotFolder  = "./Figures/GroupLevel";     % figure path as a string

% store labels for plotting
cond        = {'High Reward', 'Low Reward'};    % unconditioned stimuli labels
n.cond      = length(cond);                     % no. of UCS
paramLabel  = {'alpha', 'beta'};                % parameter labels
n.param     = length(paramLabel);               % no. of parameters

% ================== Model information ====================================
nMod        = 5;        % number of models

% specify the bounds, and the number of bins for grid search
bounds  = [0 1; % alpha
    0 50]; % beta
n.bin   = [20 30] ;

% parameter bounds [lower; upper] * parameters [b epsilon alpha(RW) beta(RW) alpha(CK) beta(CK) alpha_c beta_c]
pbounds = [0 0 0.01 0 0.01 0 0.01 0;
           1 1 1 50 1 50 1 50];
% pbounds = [0.25 0 0.01 0 0.01 0 0.01 0;
%            0.75 1 1 50 1 50 1 50];

nameModels = {'RW model' 'RW model', 'RR model', 'WSLS model',...
    'RW-CK', 'RW-CK', 'RW-CK', 'RW-CK'};

% ================== Plot colors ==========================================
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

% ================== Add paths ============================================

% add path of the current folder
tmp = fileparts(which('Step6_plotData_GroupLevel'));
addpath(tmp);

% add path to required folders in the current folder
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'HelperFunctions'));
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'LikelihoodFunctions'))
addpath(fullfile(tmp, 'SimulationFunctions'))
addpath(fullfile(tmp, 'FittingFunctions'))

% add path to the data folder
rootdir     = tmp(1:end-32);
datapath    = fullfile(rootdir, '01_Data', 'LSim_3_behavioural');

%% Section 2a: For each subject, load and analyse simulation task data

for i = 1:length(subjects)  

    [data_subj, BIC(i,:), iBEST(i), BEST(i,:), pars, NegLL(i,:)] = subjectAnalysis_v2(subjects(i), datapath, nMod, pbounds);
    
    % store subject data
    data.score(i) = data_subj.score;
    data.choice(:,i) = data_subj.choice;
    
    % store parameter values
    p(1).value(i) = pars(1,1);  % b
    p(2).value(i) = pars(2,1);  % epsilon
    p(3).value(i) = pars(3,1);  % alpha
    p(4).value(i) = pars(3,2);  % beta
    p(5).value(i) = pars(4,1);  % alpha_n
    p(6).value(i) = pars(4,2);  % beta_n
    p(7).value(i) = pars(4,3);  % alpha_c
    p(8).value(i) = pars(4,4);  % beta_c
    
end

% Store parameters in a table
alpha = p(3).value';
beta = p(4).value';
alpha_n = p(5).value';
beta_n = p(6).value';
alpha_nc = p(7).value';
beta_nc = p(8).value';
epsilon = p(2).value';
b = p(1).value';

% store in tables
t1 = table(subjects', alpha, beta, alpha_n, beta_n, alpha_nc, beta_nc, epsilon, b);
t1.Properties.VariableNames([1]) = {'Subjects'};

% Store loglikelihoods and BIC values in a table
t2 = table(subjects', NegLL);
t2.Properties.VariableNames([1]) = {'Subjects'};

t3 = table(subjects', BIC);
t3.Properties.VariableNames([1]) = {'Subjects'};
 
% Store accuracy values
t4 = table(subjects', data.score');
t4.Properties.VariableNames([1]) = {'Subjects'};
t4.Properties.VariableNames([2]) = {'Accuracy'};

% Save choice behaviour
t5 = table(subjects', data.choice');
t5.Properties.VariableNames([1]) = {'Subjects'};

% % calculate LRT
% df = length(subjects)*2;
% [h,pValue,stat,cValue] = lratiotest(LL_sum(:,3), LL_sum(:,1), df)
% t6 = table(subjects', h, pValue, stat, cValue);

% save tables
writetable(t1, sprintf("DataOutput/paramVal%i-%i.csv", min(subjects), max(subjects)));
writetable(t2, sprintf("DataOutput/NegLLVal%i-%i.csv", min(subjects), max(subjects)))
writetable(t3, sprintf("DataOutput/BICVal%i-%i.csv", min(subjects), max(subjects)))
writetable(t4, sprintf("DataOutput/accuracy%i-%i.csv", min(subjects), max(subjects)))
writetable(t5, sprintf("DataOutput/choicebehaviour%i-%i.csv", min(subjects), max(subjects)))

%% Section 2b: For each subject, load and analyse bandit task data

for i = 1:length(subjects)  

    [data_subj_band, BIC_band(i,:), iBEST_band(i), BEST_band(i,:), pars_band, NegLL_band(i,:)] = subjectAnalysis_bandit(subjects(i), datapath, nMod, pbounds);
    
    % store subject data
    data_band.score(i) = data_subj_band.score;
    data_band.choice(:,i) = data_subj_band.choice;
    
    % store parameter values
    p(1).band(i) = pars_band(1,1);  % b
    p(2).band(i) = pars_band(2,1);  % epsilon
    p(3).band(i) = pars_band(3,1);  % alpha
    p(4).band(i) = pars_band(3,2);  % beta
    p(5).band(i) = pars_band(4,1);  % alpha_n
    p(6).band(i) = pars_band(4,2);  % beta_n
    p(7).band(i) = pars_band(4,3);  % alpha_c
    p(8).band(i) = pars_band(4,4);  % beta_c
    
end

% Store parameters in a table
alpha_band = p(3).band';
beta_band = p(4).band';
alpha_n_band = p(5).band';
beta_n_band = p(6).band';
alpha_nc_band = p(7).band';
beta_nc_band = p(8).band';
epsilon_band = p(2).band';
b_band = p(1).band';

% store in tables
t7 = table(subjects', alpha_band, beta_band, alpha_n_band, beta_n_band, alpha_nc_band, beta_nc_band, epsilon_band, b_band);
t7.Properties.VariableNames([1]) = {'Subjects'};

% Store loglikelihoods and BIC values in a table
t8 = table(subjects', NegLL_band);
t8.Properties.VariableNames([1]) = {'Subjects'};

% save tables
writetable(t7, sprintf("DataOutput/paramVal-bandit%i-%i.csv", min(subjects), max(subjects)));
writetable(t8, sprintf("DataOutput/NegLLVal-bandit%i-%i.csv", min(subjects), max(subjects)))

%% Section 3a: Plot choice behaviour for simulation task
% Is the data.choice behaviour evolving over the trials to move towards the 
% p(reward|HR) threshold?

% ================== Plot p(correct data.choices)  =========================

% plot using the barScatter.m function
fh.score = figure('Name','score'); set(fh.score,'position',[100 500 250 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');

barScatter(data.score,[],[],true);   % barplot function with error bars
set(gca,'xtick',[]);
ylabel('p(HR)');
title(sprintf('subject %i-%i',subjects(1), subjects(end)), 'FontWeight','Normal');
ylim([0 1])
box off  
set(gca, 'fontsize', 14)

% save figure
if savePlots
    fh.score.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Prob_HR/prob_%i-%i.png', plotFolder,...
        min(subjects),max(subjects)))
end

% ======== Plot raw behavioural data & smoothed response function ========

% calculate HR choice mean over each trial
HR_choiceMean = nanmean(data.choice, 2);
    
% rescale the choice mean such that it ranges between 0 (LR) and 1 (HR)
if highRewAction == 2
    HR_choiceMean = (2 - HR_choiceMean)';
else
    HR_choiceMean = (HR_choiceMean - 1)';
end

% no. of trials over which to smooth
smoothingkernel = 6;

% start plot
fh.choice = figure('Name', 'Trial-by-trial choice');
box off; hold on;
set(fh.choice,'position', [500 500 700 450],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');
set(gca, 'fontsize', 14)

% set y axis limits
ylim([-0.1 1.1]);

% add data to the plot
line([0, length(HR_choiceMean)],...
    [rprob(1), rprob(1)],...
    'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
line([0, length(HR_choiceMean)],...
    [rprob(2),rprob(2)],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
plot(HR_choiceMean', ':', 'color', AZred, 'linewidth',0.5)
plot(mySmooth(HR_choiceMean, smoothingkernel,[], 'backward'),...
        '-','color', AZred,'linewidth',1); 

% add labels
title(sprintf('Mean choice, subject %i-%i',subjects(1), subjects(end)), 'FontWeight','Normal');
ylabel('p(HR)', 'fontweight','bold','fontsize',18);
xlabel('trial', 'fontweight','bold','fontsize',18);
legend({'p(reward|HR)', 'p(reward|LR)','mean choice','smoothed mean choice (HR)'},'location','southeast')
%legend boxoff

if savePlots
    fh.choice.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/choices_%i-%i.png',...
        plotFolder, min(subjects), max(subjects)))
end

%% Section 3b: Plot choice behaviour for bandit task
% Is the data.choice behaviour evolving over the trials to move towards the 
% p(reward|HR) threshold?

% ================== Plot p(correct data.choices)  =========================

% plot using the barScatter.m function
fh.score_band = figure('Name','bandit-score'); set(fh.score_band,'position',[100 500 250 400],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');
barScatter(data_band.score,[],[],true);   % barplot function with error bars
set(gca,'xtick',[]);
ylabel('p(HR)');
title(sprintf('Subject %i-%i', subjects(1), subjects(end)), 'FontWeight','Normal');
ylim([0 1])
box off  

% save figure
if savePlots
    fh.score_band.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Prob_HR/BanditTask/probBandit_%i-%i.png', plotFolder,...
        min(subjects),max(subjects)))
end

% ======== Plot raw behavioural data & smoothed response function ========

% calculate HR choice mean over each trial
HR_choiceMean = nanmean(data_band.choice, 2);
    
% rescale the choice mean such that it ranges between 0 (LR) and 1 (HR)
if highRewAction == 2
    HR_choiceMean = (2 - HR_choiceMean)';
else
    HR_choiceMean = (HR_choiceMean - 1)';
end

% no. of trials over which to smooth
smoothingkernel = 6;

% start plot
fh.choice_band = figure('Name', 'Trial-by-trial choice');
box off; hold on;
set(fh.choice_band,'position', [500 500 700 450],'paperunits','centimeters',...
    'paperposition',[0 0 6 6],'Color','w');
set(gca, 'fontsize', 12)

% set y axis limits
ylim([-0.1 1.1]);

% add data to the plot
line([0, length(HR_choiceMean)],...
    [rprob(1), rprob(1)],...
    'LineStyle', '--', 'Color', AZsky, 'linewidth',0.5);  hold on    
line([0, length(HR_choiceMean)],...
    [rprob(2),rprob(2)],...
    'LineStyle', '--', 'Color', AZcactus, 'linewidth',0.5);  hold on
plot(HR_choiceMean', ':', 'color', AZred, 'linewidth',0.5)
plot(mySmooth(HR_choiceMean, smoothingkernel,[], 'backward'),...
        '-','color', AZred,'linewidth',1); 

% add labels
title(sprintf('Mean choice, subject %i-%i',subjects(1), subjects(end)), 'FontWeight','Normal');
ylabel('p(HR stimuli)');
xlabel('trial');
legend({'p(reward|HR)', 'p(reward|LR)','mean choice','smoothed mean choice (HR)'},'location','southeast')
%legend boxoff

if savePlots
    fh.choice_band.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Choices/BanditTask/choicesBandit_%i-%i.png',...
        plotFolder, min(subjects), max(subjects)))
end

% done