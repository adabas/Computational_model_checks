function v = M3_valueUpdate(alpha, q, r, t, pt)

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

% compute value update
v = q + alpha * d;

end