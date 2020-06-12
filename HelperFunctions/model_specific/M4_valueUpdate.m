function [v, k, d] = M4_valueUpdate(alpha, alpha_c, q, k, r, t, pt)

% M3_VALUEUPDATE
% The function updates choice value based on the Rescorla Wagner rule.
%
% INPUT
%       alpha   : learning rate (alpha) parameter value
%       q       : choice value
%       r       : reward value
%       t       : trial number
%       pt      : vector containing range of partial trial values
%
% OUTPUT
%       v       : the updated choice value
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% compute prediction error
if ismember(t, pt)
    d = 0;  % no update
else
    d = r - q;
end

% update values
v = q + alpha * d;

% update choice kernel
k = (1-alpha_c) * k;
k = k + alpha_c * 1;

% % compute value update
% v = q + alpha * d;

end 