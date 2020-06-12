function [data, BIC, iBEST, BEST, pars, LL, ch] = subjectAnalysis(subID, datapath)



%% Section 1: Load data

% add path to the subjects's data folder
addpath(fullfile(datapath, sprintf('subject_0%i', subID)));

% load 'results' matrix
load(sprintf('subject_0%i_sessionII_results.mat',subID), 'results');

% load settings
load(sprintf('subject_0%i_settingsS2.mat',subID));

% load stimuli
load(sprintf('subject_0%i_sessionII_items.mat', subID), 'itemmat');

%% Section 2: Data restructuring and organisation

% store data.choices: HR coded as 1 and LR coded as 2
data.choice = cell2mat(results(2:end,9));
ch = data.choice;

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

%% Section 3: Fit models to the data

% % Find the parameter values of RW+softmax function that maximise the
% % likelihood of the data. As a sanity check, run both fmincon and grid
% % search methods
% 
% % Run fmincon optimiser
% % =====================================================================
% % set fmincon options: 
%     % - 'MaxFunEval': maximum number of function evaluations allowed.
%     % - 'Display' set to 'notify': displays ouptput if the function does
%     %                               not converge.
%     % - 'algorithm' set to 'active-set': can take large steps, which adds
%     %                               speed. The algorithm is effective on
%     %                               some problems with nonsmooth constraints.
% options=optimset('MaxFunEval', 100000, 'Display', 'notify', ...
%     'algorithm', 'active-set');
% 
% % run optimization over 10 starting points
% for count = 1:10
% 
%     % random starting points [alpha beta]
%     X0  = [rand exprnd(10)];
% 
%     % create RW+softmax function handle to input in fmincon
%     obFunc = @(x) lik_M3RescorlaWagner_v1(data.choice', data.rate.norm', x(1), x(2), idPartial);
% %         obFunc = @(x) lik_M3RescorlaWagner_v1(data.choice', data.rate.norm, X0(1), X0(2), idPartial);
% 
%     % store the lower and upper bounds of [alpha beta] parameters
%     LB = [bounds(1,1) bounds(2,1)];
%     UB = [bounds(1,2) bounds(2,2)];
% 
%     % run fmincon to check the best fitting parameter to the data
%     [xf, NegLL] = fmincon(obFunc, X0, [], [], [], [], LB, UB, [], options);
% 
%     % store starting value (Xint), fitted value (Xfit) and the negative
%     % LL values (negLL)
%     fminX.initiate(1,count) = X0(1);
%     fminX.initiate(2,count) = X0(2);
%     fminX.fit(1,count)      = xf(1);
%     fminX.fit(2,count)      = xf(2);
%     fminX.negLL(count)      = NegLL;
% 
%     % clear repeating variables from the workspace
%     clear X0 xf NegLL LB UB
% 
% end
% 
% % find global best
% [mf,i]=min(fminX.negLL(:));
% fminX.pars = fminX.fit(:,i);

%% Section 4: Save best fitting parameter and the BIC for each of the models

 [BIC, iBEST, BEST, pars, LL] = fit_all_v1(data.choice, data.rate.norm, idPartial);


end