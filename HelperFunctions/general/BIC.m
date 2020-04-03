function b = BIC(xl, cl, NegLL)

% BIC
% Function to calculate the BIC score.
%
% INPUT:
%       xl      : length of parameters
%       cl      : number of choice trials
%       NegLL   : the negative log likelihood value of the best fitting
%                 parameter
%
% OUTPUT:
%       b       : BIC value
% =========================================================================

% calculate the BIC
b = xl * log(cl) + 2*NegLL;

end