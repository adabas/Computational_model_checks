function [v, d] = M3_valueUpdate(alpha, q, r)

% M3_VALUEUPDATE
% The function updates choice value based on the Rescorla Wagner rule.
%
% INPUT
%       alpha   : learning rate (alpha) parameter value
%       q       : choice value
%       r       : reward value
%
% OUTPUT
%       v       : the updated choice value
%       d       : prediction error
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% compute prediction error
d = r - q;

% compute value update
v = q + alpha * d;

end