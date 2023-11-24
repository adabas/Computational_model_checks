%smooth Fitting data to the following models.
%       Model 1: Null model
%               parameter fixed at 0.5.
%       Model 2: Win-stay-lose-shift
%               epsilon : chooses the option probabilistically that is
%               rewarded and switches away from unrewarded
%       Model 3: Rescorla Wagner
%               alpha : learning rate
%               beta : inverse temperature function
%       Model 4: Rescorla-Wagner + choice kernel
%               alpha
%               beta
%               alpha_c : learning rate of choice kernel
%               beta_c  : inverse temperature
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

%% Section 1: Preparation
close all;
clearvars

% set seed
rng(244);

% ================== Modify ===============================================
subjects    = [11:55 57:60];
trialSeq    = 1:96;
savePlots   = 1; 
saveData    = 1;     
highRewAction = 2;      % 2 plot HR choices and 1 plot LR choices
rprob       = [0.8 0.3];
plotFolder  = './Figures/GroupLevel';

% ================== Model information ====================================
% parameter bounds [lower; upper] * parameters [b epsilon alpha(RW) beta(RW) alpha(CK) beta(CK) alpha_c beta_c]
pbounds = [0.5 0 0.05 0 0.05 0;     % parameter bounds updated to empirical data     
  0.5 1 1 25 1 25];

modNames    = {'Null', 'WSLS', 'RW', 'RW-CK', 'CK'}; % don't change the order
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
tmp = strsplit(fileparts(which('Step5_fitModel_groupLevel')), '/');

% add path to required folders in the current folder
addpath(genpath('./HelperFunctions'))
addpath('./LikelihoodFunctions')
addpath('./SimulationFunctions')
addpath('./FittingFunctions')

% add path to the data folder
rootdir     = fullfile('/', tmp{2}, tmp{3});
datapath    = fullfile(rootdir, '01_Data', 'LSim_3_behavioural');
clear tmp

% add path to VBA toolbox
if ~exist('VBA-toolbox', 'dir')
    addpath(genpath('/data/u_dabas_software/VBA-toolbox'))
end

%% Section 2: Load simulation task data and estimate model fit

% initiate structures
data    = struct();
p       = struct();

for i = 1:length(subjects)  
    
    % load and estimate model fit
    [data_subj, BIC(i,:), iBEST(i), BEST(i,:), pars, NegLL(i,:)] = subjectAnalysis(subjects(i), datapath, nMod, pbounds);
    
    % counter for storing stimuli presentation
    a = 2*i - 1; 
    
    % store subject data
    data.score(i) = data_subj.score;
    data.choice(:,i) = data_subj.choice;    % HR or LR
    data.stim(:,i) = data_subj.stimuli;     % raw choices
    data.pres(:,a:a+1) = data_subj.stimPresented;
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
t2 = table(subjects', 'VariableNames', {'Subject'});
t2 = addvars(t2, NegLL(:,1), NegLL(:,2), NegLL(:,3), NegLL(:,4), NegLL(:,5),...
    'NewVariableNames', cellfun(@(x) ['NegLL_' x], modNames, 'UniformOutput', false));

t3 = table(subjects', iBEST', 'VariableNames',  {'Subjects', 'winModel'});
t3 = addvars(t3, BIC(:,1), BIC(:,2), BIC(:,3), BIC(:,4), BIC(:,5),...
    'NewVariableNames', cellfun(@(x) ['BIC_' x], modNames, 'UniformOutput', false));
 
% Store accuracy
t4 = table(subjects', data.score', 'VariableNames', {'Subjects', 'Accuracy'});

% Save choice behaviour
t5 = table(subjects', data.stim', 'VariableNames', {'Subjects', 'Selected_stimulus'});
t6 = table(subjects', data.binRate', 'VariableNames', {'Subjects', 'Binarised_rewards'});

% Compute exceedance probability and model frequency using VBA toolbox
[posterior,out] = VBA_groupBMC(-0.5*BIC'); % BIC / AIC need to be rescaled and reversed when using the VBA toolbox (see:https://muut.com/vba-toolbox#!/vba-toolbox/questions:aicbic-calculation)
f = out.Ef;
EP = out.ep;
t7 = table(modNames', EP', f, 'VariableNames', {'Model', 'ExceedanceProbability', 'ModelFrequeny'});

% save tables
if saveData
    writetable(t1, sprintf("Data/paramVal.csv"));
    writetable(t2, sprintf("Data/NegLLVal.csv"));
    writetable(t3, sprintf("Data/BICVal.csv"));
    writetable(t4, sprintf("Data/accuracy.csv"));
    writetable(t5, sprintf("Data/choices.csv"));
    writetable(t6, sprintf("Data/rewards.csv"));
    writetable(t7, sprintf("Data/modelComparison.csv"))
end

%% Section 3: Plot smoothed choice behaviour for simulation task
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

% estimate the mean HR choice and AUC of choices for each model
m = 1;
n = 2;
for i = 1:numel(subjects)
      
    % select subject data
    stim = data.stim(:, i);
    binRate = data.binRate(:, i);
    stimPres = data.pres(:,m:n);

    % estimate the 5 models
    [~, ~, tmp.p1] = lik_M1random(stim, b(i), stimPres);
    [~, ~, tmp.p2] = lik_M2WSLS(stim, binRate, epsilon(i), stimPres);
    [~, tmp.p3, d, tmp.q3] = lik_M3RescorlaWagner(stim, binRate,...
        alpha_rw(i), beta_rw(i), stimPres);
    [~, tmp.p4] = lik_M4RWCK(stim, binRate, alpha_rwck_rw(i), beta_rwck_rw(i),...
        alpha_rwck_ck(i), beta_rwck_ck(i), stimPres);
    [~, tmp.p5] = lik_M5ChoiceKernel(stim, alpha_ck(i), beta_ck(i), stimPres);

    % calculate AUC for RW model
    auc_rw(i) = AUC(tmp.p3(:,1), trialSeq);
    
    % calculate AUC for choices
    auc_data(i) = AUC(2-data.choice(:,i), trialSeq);
    
    % replace missed trials with mean of nearest neighbours
    for j = 1:5
        tmp2 = tmp.(sprintf('p%d', j))(:,1);
        tmp2 = interp1(trialSeq(~isnan(tmp2)), tmp2(~isnan(tmp2)), trialSeq);
        PP.(sprintf('m%d', j))(i,:) = tmp2;
        
        clear tmp2
    end
    
    % store choice values for RW model
    meanQQ = nanmean(tmp.q3);
    QQ(i,:) = [mean(meanQQ(:,1:2)) mean(meanQQ(:,3:4))];
    
    m = m + 2;
    n = n + 2;
    
    clear tmp id stim stimPres binRate
end

% ---
% calculate mean and standard error of mean across subject for HR choice
for j = 1:5
    PP.(sprintf('mean%d', j)) = nanmean(PP.(sprintf('m%d', j)), 1);
    PP.(sprintf('smean%d',j)) = mySmooth(PP.(sprintf('mean%d', j)),5,[],'center');
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
semChoice = nanstd(choice,[],1)/sqrt(size(choice,1));

% smooth the mean
smoothMeanChoice = mySmooth(meanChoice,5,[],'center');

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
    t10 = table(meanChoice', semChoice', smoothMeanChoice', PP.smean1', PP.sem1', PP.smean2', PP.sem2',...
        PP.smean3', PP.sem3', PP.smean4', PP.sem4', PP.smean5', PP.sem5');
    t10.Properties.VariableNames = {'meanChoice', 'semChoice', 'smoothMeanChoice', 'nullMeanEst', 'nullSEMEst', 'wslsMeanEst', 'wslsSEMEst',...
        'rwMeanEst', 'rwSEMEst', 'rwckMeanEst', 'rwckSEMEst', 'ckMeanEst', 'ckSEMEst'};
    
    % store RW computed Q values
    t11 = table(subjects', QQ(:,1), QQ(:,2));
    t11.Properties.VariableNames = {'subId', 'Q_HR', 'Q_LR'};
    
    % store AUC_g computed using the RW computed choice probabilities
    t12 = table(subjects', auc_rw', auc_data');
    t12.Properties.VariableNames = {'subId', 'AUC_rw', 'AUC_data'};
    
    % save
    writetable(t10, sprintf("Data/realEstimatedChoices_centeredSmooth.csv"));   
    writetable(t11, sprintf("Data/LSim_Q.csv"));   
    writetable(t12, sprintf("Data/LSim_AUCg.csv"));
    
end

%% Section 4: Plot subject choice behaviour for simulation task

fh.subjAvg = figure('Name','Subject average choices and model estimates'); 
set(fh.subjAvg,'position', [500 500 700 400],'paperunits','centimeters',...
        'paperposition',[0 0 5 5],'Color','w');
set(gca, 'fontsize', 12)
hold on
box off; hold on;
ylim([0.3 1]);

% plot
pl = boundedline(trialSeq, meanChoice, semChoice, 'alpha','cmap',AZblack); hold on
hl3 = plot(PP.mean3, '-', 'color', [0 0.4470 0.7410], 'linewidth', 1.5); 

% add labels and legend
ylabel('choice');
xlabel('trial');
Lgnd = legend([pl hl3], 'real choice', modNames{3},'location', 'best');
Lgnd.Box = 'off';
title('Average real and estimated choices');

% save plot
if savePlots
    fh.subjAvg.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Model_RW.png', plotFolder))
end

% done