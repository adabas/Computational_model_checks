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
subjects    = [11:55 57:60]; %[11:18 20:21 23:55 57:60];  % specify subject IDs
trialSeq    = 1:96;
savePlots   = false;    % true will save plots in plotFolder directory
saveData    = false;     
highRewAction = 2;      % set 2 and 1 to plot HR choices and LR choices, respectively.
rprob       = [0.8 0.3];
plotFolder  = './Figures/GroupLevel';     % figure path as a string

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
pbounds = [0.5 0 0 0 0 0;     % parameter bounds updated to empirical data     
  0.5 1 1 400 1 250];

% ================== Models ===============================================
modNames    = {'RR', 'WSLS', 'RW', 'RW-CK', 'CK'}; % don't change the order
nMod        = numel(modNames);

% ================== Plot colors ==========================================
AZred   = [171,5,32]/256;
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;
AZblack = [0, 0, 0];
AZpurple = [0.4940 0.1840 0.5560];
AZgreen = [0.4660 0.6740 0.1880];

% ================== Add paths ============================================

% add path of the current folder
tmp = fileparts(which('Step6_plotData_GroupLevel_v2'));
addpath(tmp);

% add path to required folders in the current folder
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'HelperFunctions'));
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'LikelihoodFunctions'))
addpath(fullfile(tmp, 'SimulationFunctions'))
addpath(fullfile(tmp, 'FittingFunctions'))

% add path to the data folder
rootdir     = tmp(1:end-37);
datapath    = fullfile(rootdir, '01_Data', 'LSim_3_behavioural');

clear tmp

%% Section 2: For each subject, load and analyse simulation task data

% start counter
a = 1;
b = 2;

for i = 1:length(subjects)  

    [data_subj, BIC(i,:), iBEST(i), BEST(i,:), pars, NegLL(i,:)] = subjectAnalysis_v2(subjects(i), datapath, nMod, pbounds);
    
    % for missed trials, re  
    if isnan(data_subj.choice(end)); data_subj.stimuli(end+1) = 0; end
    
    % store subject data
    data.score(i) = data_subj.score;
    data.choice(:,i) = data_subj.choice;
    data.stim(:,i) = data_subj.stimuli;
    data.pres(:,a:b) = data_subj.stimPresented;
    data.binRate(:,i) = data_subj.rate.binary;
    
    % store parameter values
    p(1).value(i) = pars(1,1);  % b
    p(2).value(i) = pars(2,1);  % epsilon
    p(3).value(i) = pars(3,1);  % alpha
    p(4).value(i) = pars(3,2);  % beta
    p(5).value(i) = pars(4,1);  % alpha_n
    p(6).value(i) = pars(4,2);  % beta_n
    p(7).value(i) = pars(4,3);  % alpha_c
    p(8).value(i) = pars(4,4);  % beta_c
    p(9).value(i) = pars(5,1);  % alpha_ck
    p(10).value(i) = pars(5,2); % beta_ck
    
    a = a + 2;
    b = b + 2;
    
end

% Store parameters in a table
alpha_rw = p(3).value';
beta_rw = p(4).value';
alpha_rwck_rw = p(5).value';
beta_rwck_rw = p(6).value';
alpha_rwck_ck = p(7).value';
beta_rwck_ck = p(8).value';
alpha_ck = p(9).value';
beta_ck = p(10).value';
epsilon = p(2).value';
b = p(1).value';

% store parameter values
t1 = table(subjects', alpha_rw, beta_rw, alpha_rwck_rw, beta_rwck_rw,...
    alpha_rwck_ck, beta_rwck_ck, alpha_ck, beta_ck, epsilon, b);
t1.Properties.VariableNames(1) = {'Subjects'};

% Store negative loglikelihoods and BIC values
t2 = table(subjects', NegLL);
t2.Properties.VariableNames(1) = {'Subjects'};

t3 = table(subjects', iBEST', BIC);
t3.Properties.VariableNames(1:2) = {'Subjects', 'winModel'};
 
% Store accuracy
t4 = table(subjects', data.score');
t4.Properties.VariableNames(1) = {'Subjects'};
t4.Properties.VariableNames(2) = {'Accuracy'};

% Save choice behaviour
t5 = table(subjects', data.choice');
t5.Properties.VariableNames(1) = {'Subjects'};

% % calculate LRT
% df = length(subjects)*2;
% [h,pValue,stat,cValue] = lratiotest(LL_sum(:,3), LL_sum(:,1), df)
% t6 = table(subjects', h, pValue, stat, cValue);

% save tables
if saveData
    writetable(t1, sprintf("DataOutput/paramVal.csv"));
    writetable(t2, sprintf("DataOutput/NegLLVal.csv"));
    writetable(t3, sprintf("DataOutput/BICVal.csv"));
    writetable(t4, sprintf("DataOutput/accuracy.csv"));
    writetable(t5, sprintf("DataOutput/choices.csv"));
end


%% Section 3: Plot choice behaviour for simulation task
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
    
    if ~exist(plotFolder, 'dir'); mkdir(plotFolder); end
    
    fh.score.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/probHR.png', plotFolder))
end

% ======== Plot raw behavioural data & smoothed response function ========

% estimate the mean HR choice for each model
m = 1;
n = 2;
for i = 1:numel(subjects)
    
    % remove missed trials
    id = ~isnan(data.binRate(:,i));
    stim = data.stim(id, i);
    binRate = data.binRate(id, i);
    stimPres = data.pres(id,m:n);

    % estimate the 5 models
    [~, ~, tmp.p1] = lik_M1random_v2(stim, b(i), stimPres);
    [~, ~, tmp.p2] = lik_M2WSLS_v2(stim, binRate, epsilon(i), stimPres);
    [~, tmp.p3, d] = lik_M3RescorlaWagner_v2(stim, binRate,...
        alpha_rw(i), beta_rw(i), [], stimPres);
    [~, tmp.p4] = lik_M4RWCK_v2(stim, binRate, alpha_rwck_rw(i), beta_rwck_rw(i),...
        alpha_rwck_ck(i), beta_rwck_ck(i), [], stimPres);
    [~, tmp.p5] = lik_M5ChoiceKernel_v2(stim, alpha_ck(i), beta_ck(i), stimPres);

    % replace missed trials with mean of nearest neighbours
    for j = 1:5
        tmp2 = nan(size(data.stim(:,i)));
        tmp2(id) = tmp.(sprintf('p%d', j))(:,1);
        tmp2 = interp1(trialSeq(~isnan(tmp2)), tmp2(~isnan(tmp2)), trialSeq);
        PP.(sprintf('m%d', j))(i,:) = tmp2;
        
        clear tmp2
    end
    
    m = m + 2;
    n = n + 2;
    
    clear tmp id stim stimPres binRate
end

% ---
% calculate mean and standard error of mean across subject for HR choice
for j = 1:5
    PP.(sprintf('mean%d', j)) = nanmean(PP.(sprintf('m%d', j)), 1);
    PP.(sprintf('smean%d',j)) = mySmooth(PP.(sprintf('mean%d', j)),6,[],'backward');
    PP.(sprintf('std%d', j)) = nanstd( PP.(sprintf('m%d', j)) );
    PP.(sprintf('sem%d', j)) = nanstd( PP.(sprintf('m%d', j)) ) / sqrt( numel(subjects) );
end

% calculate the mean of real choices
if highRewAction == 2
    choice = (2 - data.choice)';
else
    choice = (data.choice - 1)';
end
meanChoice = nanmean(choice, 1);

% smooth the mean
smoothMeanChoice = mySmooth(meanChoice,6,[],'backward');

% ---
% plot the mean with shaded standard error of mean for each model
% estimation

% initialise plot
fh.est = figure('Name','Model estimates'); 
set(fh.est,'position', [500 500 700 400],'paperunits','centimeters',...
        'paperposition',[0 0 5 5],'Color','w');
set(gca, 'fontsize', 12)
hold on
box off; hold on;
ylim([0.3 1]);

% plot
pl = plot(smoothMeanChoice, '-','color', AZblack,'linewidth',1.5); hold on
hl1 = boundedline(trialSeq,PP.smean1, PP.sem1,'alpha','cmap',AZsky); hold on
hl2 = boundedline(trialSeq,PP.smean2, PP.sem2,'alpha','cmap',[0.9290 0.6940 0.1250]); hold on
hl3 = boundedline(trialSeq,PP.smean3, PP.sem3,'alpha','cmap',[0 0.4470 0.7410]); hold on
hl4 = boundedline(trialSeq,PP.smean4, PP.sem4,'alpha','cmap',AZred); hold on
hl5 = boundedline(trialSeq,PP.smean5, PP.sem5,'alpha','cmap',AZgreen); hold on

% add labels and legend
ylabel('choice');
xlabel('trial');
Lgnd = legend([pl hl1 hl2 hl3 hl4 hl5], 'real choice', modNames{1}, modNames{2}, modNames{3},...
    modNames{4}, modNames{5},'location', 'best');
Lgnd.Box = 'off';
title('Smoothed mean real and estimated choices');

% save plot
if savePlots
    fh.est.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Model_est.png', plotFolder))
end

% save data
if saveData
    
    % convert estimated choices and estimated standard error of mean to a
    % table
    t10 = table(meanChoice', smoothMeanChoice', PP.smean1', PP.sem1', PP.smean2', PP.sem2',...
        PP.smean3', PP.sem3', PP.smean4', PP.sem4', PP.smean5', PP.sem5');
    t10.Properties.VariableNames = {'meanChoice', 'smoothMeanChoice', 'nullMeanEst', 'nullSEMEst', 'wslsMeanEst', 'wslsSEMEst',...
        'rwMeanEst', 'rwSEMEst', 'rwckMeanEst', 'rwckSEMEst', 'ckMeanEst', 'ckSEMEst'};
    
    % save
    writetable(t10, sprintf("DataOutput/realEstimatedChoices.csv"));   
    
end
% done