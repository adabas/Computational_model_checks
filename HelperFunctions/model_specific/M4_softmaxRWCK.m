function p = M4_softmaxRWCK(q, k, beta, beta_c)

% M3_SOFTMAXFUNCTION
% Compute choice probabilities for value update and choice kernel.
%
% INPUT: 
%       q       : a vector containing choice value
%       k       : a vector containing choice kernel
%       beta    : inverse temperature parameter
%       beta_c  : choice kernel temperature parameter
%
% OUTPUT:
%       p       : vector containing choice probability
%
% Aroma Dabas [dabas@cbs.mpg.de]
% April 2020
% =========================================================================

V = (beta * q) + (beta_c * k);
p = exp(V) / sum(exp(V));


end