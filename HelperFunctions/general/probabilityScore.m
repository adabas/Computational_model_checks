function [choice, score, tr_range] = probabilityScore(tr, tr_l, sim, m, trialSelect)


% Calculate the probability of selecting high reward stimuli

if strcmp(trialSelect, 'First') % first set of trials
        
    % find trial bounds
    tr_range = [1 tr];

    % save the stimuli selection for the trials
    choice = sim(m).a(1:tr,:);

    % convert 2 to 1 and 1 to 0
    choice = choice -1;

    % calculate score
    score = (nansum(choice == 1, 1))/tr;

elseif strcmp(trialSelect, 'Middle') % middle set of trials

    % find trial bounds
    tr_range = [tr_l/2-(tr/2) tr_l/2+(tr/2)-1];

    % save the stimuli selection for the trials
    choice = sim(m).a(tr_range(1):tr_range(2),:);

    % convert 2 to 1 and 1 to 0
    choice = choice -1;

    % calculate score
    score = (nansum(choice == 1, 1))/tr;

elseif strcmp(trialSelect, 'Last') % last set of trials

    % find trial bounds
    tr_range = [(tr_l-tr)+1 tr_l];

    % save the stimuli selection for the trials
    choice = sim(m).a(tr_range(1):tr_range(2),:);

    % convert 2 to 1 and 1 to 0
    choice = choice -1;

    % calculate score
    score = (nansum(choice == 1, 1))/tr;

end

end