function s = stimuliPresentation(t)

% STIMULIPRESENTATION Function to create presentation of stimuli.
%
% INPUT:
%   t   : number of trials
%
% Output:
%   s   : t X 2 matrix of trial specific stimuli presentation
%
% Aroma Dabas
% October 2022
% -------------------------------------------------------------------------

% Create order of the stimuli pairs without repeating the same presentation
% twice
nconditions = 4; 
nrepeats = 24; %6;
n = nconditions * nrepeats;

combInfo = {'one', 'two', 'three', 'four'};

order = cell(1,n);
order(1,1:nconditions) = combInfo(randperm(nconditions)); % include the first set of combnations
for k = 1:nrepeats-1
    m = k * nconditions;
    r = combInfo(randperm(nconditions)); %randperm(nconditions);
    if strcmp(r(1), order(m)) %( r(1) == order(m) )
        r = fliplr(r);
    end
   order(1,m+1:m+nconditions) = r;
end

% ---

% intialise stimuli matrix
s = NaN(t, 2);

% Create combnations
comb.one = [1 3];
comb.two = [1 4];
comb.three = [2 3];
comb.four = [2 4];

% stimuli presentation
for i = 1:t
   r = randi([1,2],1);
    
   s(i,1) = comb.(sprintf('%s', order{i}))(r);
   if s(i,1) == comb.(sprintf('%s', order{i}))(1)
       s(i,2) = comb.(sprintf('%s', order{i}))(2);
   else
       s(i,2) = comb.(sprintf('%s', order{i}))(1);
   end
    
end

end