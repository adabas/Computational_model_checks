function [q, k, d] = M4_valueUpdate(alpha, alpha_c, q, k, r, t, pt)

% M3_VALUEUPDATE
% The function updates choice value based on the Rescorla Wagner rule.
%
% INPUT
%       alpha   : learning rate parameter value
%       alpha_c : choice kerne parameter value
%       q       : choice value
%       k       : kernel
%       r       : reward value
%       t       : trial number
%       pt      : vector containing range of partial trial values
%
% OUTPUT
%       q       : updated choice value
%       k       : updated kernel
%       d       : prediction error
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% compute whether the same or a differnt choice was selected
if ismember(t, pt)
    d = 0;  % no update for missed/partial
else
    d = r - q;
end

% update value
q = q + alpha * d;

% update choice kernel
k = (1-alpha_c) * k;
k = k + alpha_c * 1;

end 