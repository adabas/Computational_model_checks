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

 % Is the reward (rand()) uniform? Perhaps it comes from a gaussian
% distribution. There are options for testing this.
% One simple check is to use rewards from a Bernoulli distribution. For
% this, we can take the values 1 (pleasant) and 0 (not pleasant) or -1
% (unpleasant).
% Another (and probably better) approach is to estimate the reward
% distributions from the empirical data, and use this distribution for the
% simulated data. Additionally, we can check if the estimated distribution
% is similar to a gaussian distribution.
 
 end 