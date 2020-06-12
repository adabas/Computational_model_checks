function p = M4_softmaxCK(q, k, beta, beta_c)

% M3_SOFTMAXFUNCTION
% Softmax function
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

V = beta * q + beta_c * k;
p = exp(V) / sum(exp(V));


end