function [data, BIC, iBEST, BEST, pars, NegLL] = subjectAnalysis_bandit(subID, datapath, nMod, pbound)

% SUBJECTANALYIS
% This function loads subject information and restructures the data. Also,
% it fits data to the models.

%% Section 1: Load data

% add path to the subjects's data folder
addpath(fullfile(datapath, sprintf('subject_00%i', subID)));

% load results
load(sprintf('subject_%i_bandit_results.mat',subID), 'results');%, 'settings');

%% Section 2: Data restructuring and organisation

% store choices
data.choice = zeros((length(results)-1),1);
for i = 1:(length(results)-1)
    if strcmp(results{i+1,6}, 'NAN')
        data.choice(i) = NaN;
    else
        data.choice(i) = results{i+1,6};
    end
end

% store the number of trials
n.trials = length(data.choice);

% calculate the probability with which the participant selected HR stimuli
data.score = (sum(data.choice == 1))/(n.trials - sum(isnan(data.choice)));

% load outcome
for i = 1:n.trials
    if strcmp(results{i+1,9}, 'NaN')
       data.outcome(i) = NaN;    % replace the partial trial events with NaN
    elseif strcmp(results{i+1,9}, '+1')
       data.outcome(i) = 1;
    elseif strcmp(results{i+1,9}, '-1')
       data.outcome(i) = 0;
    end
end
data.outcome = data.outcome';

% left vs right key selection
for i = 1:n.trials
    if strcmp(results{i+1,7}, 'Missed')
       keyChoice(i) = NaN;    % replace the partial trial events with NaN
    elseif strcmp(results{i+1,7}, 'Left-key')
       keyChoice(i) = 1;
    elseif strcmp(results{i+1,7}, 'Right-key')
       keyChoice(i) = 2;
    end
end
keyChoice(isnan(keyChoice)) = [];

keyChoice = keyChoice';

%% Section 4: Save best fitting parameter and the BIC for each of the models

[BIC, iBEST, BEST, pars, NegLL] = fit_all_v1(data.choice, data.outcome, [], nMod, pbound, keyChoice);

end