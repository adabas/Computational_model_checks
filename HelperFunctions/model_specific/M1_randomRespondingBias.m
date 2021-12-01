function p = M1_randomRespondingBias(b)

% M1_RANDOMRESPONDINGBIAS
% Function to compute choice probabilities based on the random responding
% model.
%
% INPUT:
%       b   : bias parameter for a choice
%
% OUTPUT:
%       p   : vector indicating choice probabilities
%
% Aroma Dabas [dabas@cbs.mpg.de]
% =========================================================================

p = [b 1-b];

end