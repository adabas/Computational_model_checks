% The script simulates a data set and estimates the parameters using both
% the grid search approach and the fmincon optimisers. This is a sanity
% check to see the difference between the two approaches.
%
% For a model with two free parameters, the grid search approach should be
% as good as the fmincon procedure for estimating the parameters. The
% difference should, however, arise when we go for models with more than
% two free parameters.
%
% Currently, it checks for one simulated data. In future, this can be
% implemented for a number of simulated data to get a more generalisable
% effect.
% 
% Aroma Dabas [dabas@cbs.mpg.de]
% Max Planck Institute for Human Cognitive and Brain Sciences
% 
%   01-2020         version 1
% =========================================================================

%% Section 1: Preparations

clearvars; close all

% set seed
rng(234, 'twister');

% ================== Modify ===============================================
% experiment parameters
T       = 100;      % number of trials
rbounds = [0 1];    % bounds of the mean reward
nRep    = 40;   % 100    % number of simulation repetitions
rprob   = [0.8 0.35]; % reward probability for [HR LR] stimuli
Npt     = 0;       % number of partial trials

% store information about the parameter
paramlabel  = {'alpha', 'beta', 'alpha_c', 'beta_c'};
n.param1    = 2;    % number of parameters in model RW
n.param2    = 4;    % number of parameters in model RW-CK

% state the bounds and bins of each parameter
bounds  = [0 1; % alpha
    1 7;    % beta
    0 1;    % alpha_c
    1 10];  % beta_c

n.bin1   = [20 30];
n.bin2   = [20 30 20 30];

% parameter setting for simulating the data for RW+softmax function
realalpha   = .2;
realbeta    = 10;
realalpha_c = 0.1;
realbeta_c  = 4;

% ================== Add paths ============================================

addpath('./SimulationFunctions')
addpath('./AnalysisFunctions')
addpath(genpath('./HelperFunctions'))
addpath('./FittingFunctions')
addpath('./LikelihoodFunctions')

% ================== Plot settings ========================================

savePlots   = 0;    % set as 1 if you want to save plots, otherwise 0
plotFolder  = './Figures/ModelSimulation/';

%% Section 2a: Simulate one data set for simple RW model

[choice.m1, reward.m1, pt.m1] = simulate_M3RescorlaWagner_v1(T, realalpha, realbeta,...
                            rprob, rbounds, Npt);

%% Section 2b: Simulate one data set for RW-CK model

[choice.m2, reward.m2, pt.m2] = simulate_M4RWCK_v1(T, realalpha, realbeta,...
    realalpha_c, realbeta_c, rprob, rbounds, Npt);
                        
%% Section 3a: fmincon optimiser for simple RW model

% set fmincon settings
options = optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% run optimization over 10 starting points
for count = 1:10
    
    % random starting points [alpha beta]
    X0  = [rand exprnd(10)];

    % create RW+softmax function handle to input in fmincon
    obFunc = @(x) lik_M3RescorlaWagner_v1(choice.m1, reward.m1, x(1), x(2), pt.m1);

    % store the lower and upper bounds of [alpha beta] parameters
    LB = [bounds(1,1) bounds(2,1)];
    UB = [bounds(1,2) bounds(2,2)];

    % run fmincon to check the best fitting parameter to the data
    [xf, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

    % store starting value (Xint), fitted value (Xfit) and the negative
    % LL values (negLL)
    fminX.initiate1(1,count) = X0(1);
    fminX.initiate1(2,count) = X0(2);
    fminX.fit1(1,count)      = xf(1);
    fminX.fit1(2,count)      = xf(2);
    fminX.negLL1(count)      = NegLL;

    % clear repeating variables from the workspace
    clear X0 xf NegLL LB UB

end

% find global best
[mf,i]=min(fminX.negLL1(:));
fminX.pars1 = fminX.fit1(:,i);

%% Section 3b: fmincon optimiser for RW-CK model

% set fmincon settings
options = optimset('MaxFunEval', 100000, 'Display', 'notify', ...
    'algorithm', 'active-set');

% run optimization over 10 starting points
for count = 1:10

   % create function capturing the RW with softmax function
    obFunc = @(x) lik_M4RWCK_v1(choice.m2, reward.m2, x(1), x(2), x(3), x(4), pt.m2);

    % create vector storing random alpha and beta starting values [alpha beta alpha_c beta_c]
    X0 = [rand exprnd(4) rand 0.5+exprnd(1)];

    % store the lower and upper bounds of [alpha beta] parameters
    LB = bounds(:,1); %[bounds(1,1) bounds(2,1)];
    UB = bounds(:,2); %[bounds(1,2) bounds(2,2)];

    % estimate the parameter fit and the neg log
    [xf, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);

    % store starting value (Xint), fitted value (Xfit) and the negative
    % LL values (negLL)
    for i = 1:length(bounds)
        fminX.initiate2(i,count) = X0(i);
        fminX.fit2(i,count)      = xf(i);
    end
    fminX.negLL2(count)      = NegLL;

    % clear repeating variables from the workspace
    clear X0 xf NegLL LB UB

end

% find global best
[mf,i]=min(fminX.negLL2(:));
fminX.pars2 = fminX.fit2(:,i);

%% Section 4a: grid search for simple RW model

% create bins for both the parameters
for iParam = 1:n.param1
    range = linspace(bounds(iParam,1),bounds(iParam,2),n.bin1(iParam)+1);
    p1{iParam} = range(2:end); % stay just off the zero bounds
end

% select each bin and determine the likelihood of the data given the
% parameters in that bin
params = nan(1,n.param1);
for t = 1:(n.bin1(1))
    params(1) = p1{1}(t);    % alpha value
    for tt = 1:(n.bin1(2))
        params(2) = p1{2}(tt);   % beta value

        % determine the log likelihood value for the parameter
        llh1(t,tt) = -lik_M3RescorlaWagner_v1(choice.m1,...
            reward.m1, params(1), params(2), pt.m1);
    end
end

%% Section 4b: grid search for Rw-CK model

% create bins for both the parameters
for iParam = 1:n.param2
    range = linspace(bounds(iParam,1),bounds(iParam,2),n.bin2(iParam)+1);
    p2{iParam} = range(2:end); % stay just off the zero bounds
end

% select each bin and determine the likelihood of the data given the
% parameters in that bin
params = nan(1,n.param2);
for t = 1:(n.bin2(1))
    params(1) = p2{1}(t);    % alpha value
    for tt = 1:(n.bin2(2))
        params(2) = p2{2}(tt);   % beta value
        for ttt = 1:(n.bin2(3))
            params(3) = p2{3}(ttt);   % alpha_c value
            for tttt = 1:(n.bin2(4))
                params(4) = p2{4}(tttt);   % beta_c value

                % determine the log likelihood value for the parameter
                llh2(t,tt,ttt,tttt) = -lik_M4RWCK_v1(choice.m2,...
                    reward.m2, params(1), params(2), params(3), params(4),pt.m2);

            end
        end
    end
end

%% Section 5a: Grid plot for RW model
% Plot parameters used for simulation with parameters estimated by grid search and fmincon optimiser

% create figure
fh = figure;

% find minimum and maximum values from the grid search
mi      = min(llh1(:));
[ma,i]  = max(llh1(:));

% create matrices containing bins of alpha and beta parameters
x = repmat(1:length(p1{1,1}), length(p1{1,2}), 1)';   % repeating 
y = repmat(1:length(p1{1,2}), length(p1{1,1}), 1);

% plot grid surface
imagesc(p1{1,1}(1:end), p1{1,2}(1:end), llh1',[mi,ma])
colorbar

% plot the parameters
hold on
plot(realalpha, realbeta, 'xr')     % simulation parameters
plot(p1{1,1}(x(i)), p1{1,2}(y(i)), 'ok')  % grid search
plot(fminX.pars1(1), fminX.pars1(2),'*k')  % fmincon optimisation

% add labels and edit settings
xlabel('alpha')
ylabel('beta')
title('RW model');
set(gca,'fontsize', 14)
dim = [0.65 0 0.5 0.05];
annotation('textbox',dim, 'String', '*: fmincon; o: grid search',...
    'EdgeColor', 'none', 'FontSize', 11);

% save plot
if savePlots
    fh.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Grid-vs-fmincon/Grid_RW_model.png', plotFolder))
end

%% Section 5b: Grid subplots for RW-CK model
% Plot parameters used for simulation with parameters estimated by grid search and fmincon optimiser

% create figure
fh2 = figure;

llh2 = squeeze(max(max(llh2,[],4),[],3));

% find minimum and maximum values from the grid search
mi      = min(llh2(:));
[ma,i]  = max(llh2(:));

% create matrices containing bins of alpha and beta parameters
x = repmat(1:length(p2{1,1}), length(p2{1,2}), 1)';   % repeating 
y = repmat(1:length(p2{1,2}), length(p2{1,1}), 1);

% plot grid surface
imagesc(p2{1,1}(1:end), p2{1,2}(1:end), permute(llh2, [1 2]), [mi,ma])
colorbar

% plot the parameters
hold on
plot(realalpha, realbeta, 'xr')     % simulation parameters
plot(p2{1,1}(x(i)), p2{1,2}(y(i)), 'ok')  % grid search
plot(fminX.pars2(1), fminX.pars2(2),'*k')  % fmincon optimisation

% add labels and edit settings
xlabel('alpha')
ylabel('beta')
title('RW-CK model');
set(gca,'fontsize', 14)
dim = [0.65 0 0.5 0.05];
annotation('textbox',dim, 'String', '*: fmincon; o: grid search',...
    'EdgeColor', 'none', 'FontSize', 11);

% save plot
if savePlots
    fh2.PaperPositionMode = 'auto';
    saveas(gcf, sprintf('%s/Grid-vs-fmincon/Grid_RW-CK_model.png', plotFolder))
end

% done