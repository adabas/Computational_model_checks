function k = M5_CKUpdate(alpha_c, k)

% M5_CKUPDATE
% The function updates choice kernel.
%
% INPUT
%       alpha   : learning rate parameter value
%       r       : reward value
%       t       : trial number
%       pt      : vector containing range of partial trial values
%
% OUTPUT
%       k       : the updated choice kernel
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

% update choice kernel
k = (1-alpha_c) * k;
k = k + alpha_c * 1;

end 