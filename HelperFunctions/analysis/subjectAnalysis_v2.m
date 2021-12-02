function [data, BIC, iBEST, BEST, pars, NegLL] = subjectAnalysis_v2(subID, datapath, nMod, pbound)

% SUBJECTANALYIS_v2
% This function loads and restructures the subject data. The restructured
% data is fed into a function that determines the fit of the model to the
% data.
%
% INPUTS:
%       subID       : subject id as double
%       datapath    : root directory of the subject folders as char var.
%       nMod        : number of models to test
%       pbound      : lower and upper parameter bounds for each model
%
% OUTPUT:

%% Section 1: Load data

% add path to the subjects's data folder
addpath(fullfile(datapath, sprintf('subject_00%i', subID)));

% load results matrix
load(sprintf('subject_%i_P3_results.mat',subID), 'results');

% load settings
load(sprintf('subject_%i_settingsP1.mat',subID));

% load stimuli
load(sprintf('subject_%i_sessionII_items.mat', subID), 'itemmat');

%% Section 2: Data restructuring and organisation

% store data.choices: HR coded as 1 and LR coded as 2
data.choice = cell2mat(results(2:end,9));

% store the number of trials
n.trials = settings.design.ntrialblock * settings.design.nblock;

% calculate the probability with which the participant selected HR stimuli
data.score = (sum(data.choice == 1))/(n.trials - sum(isnan(data.choice)));%(sum(data.choice == 1))/(n.trials);

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

% binary rewards? 0 (neutral-negative) <= 0.5; 1 (positive) > 0.5
data.rate.binary = data.rate.norm;
data.rate.binary(data.rate.norm > 0.5, 1) = 1;
data.rate.binary(data.rate.norm <= 0.5, 1) = 0;
% data.rate.binary(isnan(data.rate.norm)) = 'NaN';

% store the event types shown on each trial (Positive, Negative-Neutral & 'Partial')
data.eventType   = results(2:end, 12);
idPartial  = find(cellfun('isempty', data.eventType));  % the empty cells are partial trials
data.eventType(idPartial) = {'Partial'};

% store the corresponding condition ids of the six stimuli: HR coded as 1,
% LR coded as 2 and baseline coded as 3
for i = 1:length(itemmat)
    cond_id(i,1) = itemmat(i).cond_id;
end

% % determine key presses
% keyPress_tmp = results(2:101,8);
% data.keyPress(strcmp(keyPress_tmp, 'Left-key')) = 1;
% data.keyPress(strcmp(keyPress_tmp, 'Right-key')) = 0;
% data.keyPress = data.keyPress';

% left vs right key selection
for i = 1:n.trials
    if isempty(results{i+1,14})
       keyChoice{i} = NaN;    % replace the partial trial events with NaN
    else
       keyChoice{i} = results{i+1,8};
    end
end
keyChoiceCode = zeros(1, length(keyChoice));
keyChoiceCode(strcmp(keyChoice, 'Left-key')) = 1;
keyChoiceCode(strcmp(keyChoice, 'Right-key')) = 2;
keyChoiceCode(strcmp(keyChoice, 'Missed')) = [];

%% Section 4: Save best fitting parameter and the BIC for each of the models

[BIC, iBEST, BEST, pars, NegLL] = fit_all_v1(data.choice, data.rate.binary, idPartial, nMod, pbound, keyChoiceCode);

end