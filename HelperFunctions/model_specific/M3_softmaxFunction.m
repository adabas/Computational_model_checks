function p = M3_softmaxFunction(q, beta)

% M3_SOFTMAXFUNCTION
% Softmax function
%
% INPUT: 
%       q       : a vector containing choice value
%       beta    : inverse temperature parameter
%
% OUTPUT:
%       p       : vector containing choice probability
%
% Aroma Dabas [dabas@cbs.mpg.de]
% March 2020
% =========================================================================


ev  = exp(beta*q);
sev = sum(ev);
p = ev / sev;


end