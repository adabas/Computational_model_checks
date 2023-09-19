function [data, idPartial] = loadSubjectData(subID, datapath)

% LOADSUBJECTDATA 
% Function to load subject trial choices and rewards.
%
% INPUT:
%       subID       : binary subject id
%       datapath    : path to the folder storing all subjects' behavioural
%                     output files
%
% OUTPUT:
%       data        : structure with fields choices, stimuli presented,
%                     accuracy score, rewards, and reward type
%       idPartial   : partial trial ids
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2022
% -------------------------------------------------------------------------

% add path to the subjects's data folder
addpath(fullfile(datapath, sprintf('subject_00%i', subID)));

% load results matrix
load(sprintf('subject_%i_P3_results.mat',subID), 'results');

% load settings
load(sprintf('subject_%i_settingsP1.mat',subID));

% load stimuli
load(sprintf('subject_%i_sessionII_items.mat', subID), 'itemmat', 'stimmat');

%% Section 2: Data restructuring and organisation

% store reward selection: HR coded as 1 and LR coded as 2
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

% binarise rewards: 0 (neutral-negative) <= 0.5; 1 (positive) > 0.5
data.rate.binary = data.rate.raw; %data.rate.norm;
data.rate.binary(data.rate.raw > 0.5, 1) = 1;
data.rate.binary(data.rate.raw <= 0.5, 1) = 0;

% store the event types shown on each trial (Positive, Negative-Neutral & 'Partial')
data.feedback   = results(2:end, 12);
idPartial  = find(cellfun('isempty', data.feedback));  % the empty cells are partial trials
data.feedback(idPartial) = {'Missed'};

% store trial wise stimuli presentated and choices
for i = 1:length(itemmat)
    cond_id(i,1) = itemmat(i).ids{1};
    cond_id(i,2) = itemmat(i).cond;
end

% store trial wise stimuli selected
data.stimuliID = cell2mat(results(2:end, 7));

% loop over trials to avoid accidental stimuli ids that match cond_id
% replacements.
for i = 1:numel(data.stimuliID)
    if data.stimuliID(i) == cond_id(1,1)
        data.stimuli(i) = cond_id(1,2);
    elseif data.stimuliID(i) == cond_id(2,1)
        data.stimuli(i) = cond_id(2,2);
    elseif data.stimuliID(i) == cond_id(3,1)
        data.stimuli(i) = cond_id(3,2);
    elseif data.stimuliID(i) == cond_id(4,1)
        data.stimuli(i) = cond_id(4,2);
    elseif isnan(data.stimuliID(i))
        data.stimuli(i) = 0;
    end
end

data.stimPresented = [stimmat.presentation(1).trials; stimmat.presentation(2).trials;...
    stimmat.presentation(3).trials; stimmat.presentation(4).trials];

end