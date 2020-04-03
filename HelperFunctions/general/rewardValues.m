function r = rewardValues(rbounds)

% REWARDVALUES
% Function calculates a 1X2 vector of reward values-[highReward lowReward].
% The highReward is higher than the middle value of the reward bounds,
% whereas the lowReward is lower than the middle value of the reward
% bounds.
%
% INPUT:
%   rbounds : a 1X2 vector containing the lower and higher bounds of
%             potential reward values [lowerBound upperBound]
%
% OUTPUT:
%   r       : a 1X2 vector containing the randomly selected higher and
%             lower reward values
% =========================================================================

% middle value of the reward bounds
m = (rbounds(2) - rbounds(1)) / 2;

% randomly create vector of reward values
r = [abs(rbounds(2)-m*rand()) abs(rbounds(1)-m*rand())];

end