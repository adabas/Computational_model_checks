function a = choose(p)
%CHOOSE   This function generates the action
%
%   Input
%       p: choice probability in favor of action a = 0
%
%   Output
%       a: chosen action

% Generate action from Bernoulli distribution
a = binornd(1, p)+1;

end