function a = choose(p)
%CHOOSE   This function generates the action
%
%   Input
%       p: choice probability in favor of action a = 0
%
%   Output
%       a: chosen action


% Cumsum: cummulative sum of elements
% Rand: returns a random integer between 0 and 1
% eps:  returns the distance from 1.0 to the next largest double-precision number
%a = max(find([-eps cumsum(p)] < rand));
% Question: Why don't you use a Bernoulli distribution to choose the action? For now I added +1 but you can consider
% using A \in {0,1} instead of {1,2}
% Answer: This is a brilliant alternative! I feel the +1 does a good job.
% With [0 1], I will have to update the script to avoid using A as an index
% for selecting the reward probability (Matlab doesn't take 0 as an index
% unlike python), and not sure if that is worth it.

% Generate action from Bernoulli distribution
a = binornd(1, p)+1;

end