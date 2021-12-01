function p = M5_softmaxCK(k, beta_c)

% M5_SOFTMAXFUNCTION
% Compute choice probabilities for choice kernel.
%
% INPUT: 
%       k       : a vector containing choice kernel
%       beta_c  : choice kernel temperature parameter
%
% OUTPUT:
%       p       : vector containing choice probability
%
% Aroma Dabas [dabas@cbs.mpg.de]
% October 2020
% =========================================================================

V = beta_c * k;
p = exp(V) / sum(exp(V));


end