function auc_g = AUC(m, t)

% Compute AUC using the ground truth. We use the formula as discussed
% in Pruessner et al., 2003.
%
%   INPUT:
%       m   : 1 X n vector measurement
%       t   : 1 X n measurement time points
%
%   OUTPUT:
%       auc_g : AUC value with respective to the ground truth, here 0.
% 
% Aroma Dabas
% March 2023
% -------------------------------------------------------------------------

% remove any missed trials
imiss = isnan(m);
m(imiss) = [];
t(imiss) = [];

% ---
% distance between the measurements
d = t(2:end) - t(1:end-1);

% ---
for ii = 1:numel(m)-1
    a(ii) = (m(ii+1) + m(ii))*d(ii)/2;
end

auc_g = sum(a);

end

