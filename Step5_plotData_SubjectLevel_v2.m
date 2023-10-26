% Fitting the models to each dataset.
% 
% Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 12/03/19
% Latest update on 28/11/2022
% =========================================================================
%% Section 1: Preparation

close all;
%clearvars

rng(244, 'twister');    % set seed

% ================== Modify ===============================================
subjects    = 11:60;       % specify subject ID
savePlots   = true;     % true will save plots in plotFolder directory
plotFolder  = "./Figures/SubjectLevel";
trialSeq       = 1:96;     % for plotting

% parameter bounds [lower; upper] * parameters [b epsilon alpha(RW) beta(RW) alpha(CK) beta(CK) alpha_c beta_c]
pbounds = [0.5 0 0.05 0 0.05 0;     % parameter bounds updated to empirical data     
  0.5 1 1 25 1 25];

% ================== Add paths ============================================
% add path of the current folder
tmp = fileparts(which('Step5_plotData_SubjectLevel'));
addpath(tmp);

% add path to required folders in the current folder
addpath(genpath(fullfile(tmp, 'HelperFunctions')))
addpath(fullfile(tmp, 'LikelihoodFunctions'))
addpath(fullfile(tmp, 'SimulationFunctions'))
addpath(fullfile(tmp, 'FittingFunctions'))
addpath(genpath(fullfile(tmp, 'DataOutput')))

% add path to the data folder
tmp         = fileparts(pwd);
rootdir     = tmp(1:end-numel('03_Scripts/'));
datapath    = fullfile(rootdir, '01_Data', 'LSim_3_behavioural');

% ================== Models ===============================================
modNames    = {'RR', 'WSLS', 'RW', 'RW-CK', 'CK'}; % don't change the order
nMod        = numel(modNames);

% ================== Plot colors ==========================================
AZred   = [0.6350 0.0780 0.1840];
AZblack = [0, 0, 0];
AZcactus = [92, 135, 39]/256;
AZsky   = [132, 210, 226]/256;

clear tmp

%% BIG LOOP

for isub = subjects
    %% Section 2: Load subject data

    data = loadSubjectData(isub, datapath);
    
    % remove missed trials
    if any(data.stimuli == 0)
        [id,~] = find(data.stimuli == 0);
        data.stimuli(id,:) = [];
        data.rate.binary(id,:) = [];
        data.stimPresented(id,:) = [];
    end

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
    title(sprintf('Subject %i', isub));
    ylim([0 1])
    box off  

    % save figure
    if savePlots
        fh.score.PaperPositionMode = 'auto';
        saveas(gcf, sprintf('%s/prob_0%i.png', plotFolder, isub))
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
    
    % plot
    plot(data.choice_plot, '.r', 'MarkerSize',20, 'color', [0 0 0]);
    plot(data.choice_smooth, '-','color', AZblack,'linewidth',1.5);%[0.8500, 0.3250, 0.0980],'linewidth',1);

    % add labels
    title(sprintf('subject %i',isub));
    ylabel('choice');
    xlabel('trial');
    legend({'choice','smoothed choice'},'location','best')
    legend boxoff

    if savePlots
        fh.choice.PaperPositionMode = 'auto';
        %print -depsc2 finalPlot1.eps
        saveas(gcf, sprintf('%s/choices_0%i.png', plotFolder, isub))
    end

    %% Section 4: Model fit

    [BIC, iBEST, BEST, pars, NegLL] = fit_all_models(data.stimuli, data.rate.binary, 5, pbounds, data.stimPresented);

    %% Section 5: Plot estimated parameters
    
    % estimate the 5 models
    [~, ~, tmp.p1] = lik_M1random(data.stimuli, pars(1,1), data.stimPresented);
    [~, ~, tmp.p2] = lik_M2WSLS(data.stimuli, data.rate.binary, pars(2,1), data.stimPresented);
    [~, tmp.p3, d] = lik_M3RescorlaWagner(data.stimuli, data.rate.binary,...
        pars(3,1), pars(3,2), data.stimPresented);
    [~, tmp.p4] = lik_M4RWCK(data.stimuli, data.rate.binary, pars(4,1), pars(4,2),...
        pars(4,3), pars(4,4), data.stimPresented);
    [~, tmp.p5] = lik_M5ChoiceKernel(data.stimuli, pars(5,1), pars(5,2), data.stimPresented);
    
    % replace missed trials with mean of nearest neighbours
    for i = 1:5
        PP.(sprintf('m%d', i)) = nan(size(data.choice));
        PP.(sprintf('m%d', i))(~isnan(data.choice)) = tmp.(sprintf('p%d', i))(:,1);
        PP.(sprintf('m%d', i)) = interp1(trialSeq(~isnan(PP.(sprintf('m%d', i)))),...
            PP.(sprintf('m%d', i))(~isnan(PP.(sprintf('m%d', i)))), trialSeq);
    end
    
    % ---
    % plot the 5 models in a grid
    
    % initialise plot
    fh.est = figure('Name','Model estimates'); 
    set(fh.est,'position',[10 100 800 800],'paperunits','centimeters','Color','w');
    set(gca, 'fontsize', 12)
    
    for i = 1:5
        
        hold on
        
        % open the model's subplot
        h(i) = subplot(ceil(5/2), 2, i);
        box off; hold on;
        ylim([-0.1 1.1]);
        
        % plot
        plot(data.choice_plot, '.r', 'MarkerSize',5, 'color', [0 0 0]); hold on
        plot(data.choice_smooth, '-','color', AZblack,'linewidth',1.5); hold on
        plot(PP.(sprintf('m%d', i)),'-.','color', AZred,'linewidth',1.3);
        
        % add labels
        title(modNames{i});
        ylabel('choice');
        xlabel('trial');
        
    end
    
    % place legend
    fh.est.Position(3) = fh.est.Position(3) + 250;
    Lgnd = legend({'choice','smoothed choice', 'estimated HR choice'},...
             'location','best');
    Lgnd.Box = 'off';
    Lgnd.Position(1) = 0.65;
    Lgnd.Position(2) = 0.2;
    
    % save plot
    if savePlots
        fh.choice.PaperPositionMode = 'auto';
        saveas(gcf, sprintf('%s/Model_est_0%i.png', plotFolder, isub))
    end
    

    % Plot prediction error for RW model
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
    pe_line = plot(d,'-','color', AZsky,'MarkerEdgeColor', 'none', 'linewidth', 1.5);

    % more modifications
    set(gca, ...
        'FontName'   , 'Helvetica',...
        'TickDir'   , 'out' ,...
        'TickLength', [.01 .01],...
        'XColor'    , [.3 .3 .3],...
        'YColor'    , [.3 .3 .3],...
        'YTick'     , -0.8:0.2:0.8,...
        'XTick'     , 1:11:100,...
        'YGrid'     , 'off');

    % add labels
    ylabel('prediction error', 'FontSize', 15);
    xlabel('trial', 'FontSize', 15);

    if savePlots
        fh.delta.PaperPositionMode = 'auto';
        saveas(gcf, sprintf('%s/PredictionError_0%i.png', plotFolder, isub))
    end
    
    %% cleanup
    close all
    clear fh PP* d Q data
end
% done