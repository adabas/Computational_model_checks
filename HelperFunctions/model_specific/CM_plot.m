function CM_plot(T, mu, varargin)

% CM_plot
% Function to compute and plot a confusion matrix (CM) for the set models
% (random responding, WSLS, RW). The CM is based on the log likelihood
% values.
%
%   INPUT
%       T       : number of trials
%       mu      : the probabilistic reward structure
%
%  VARARGIN
%       pbounds : parameter values bounds
%       nMod    : number of models
%       nRep    : number of iterations
%       Npt     : number of partial trials
%       rbounds : reward bounds
%       name    : plot name
%
%  OUTPUT
%       Confusion matrix
%
% Aroma Dabas [dabas@cbs.mpg.de]
% June 2020
% =========================================================================

%% Parse input arguments

% default values
pbounds = [0.4 0.6 0.6 0.5;     % LB
    0.6 1 1 2];                 % UP
nMod    = 3;
nRep    = 100;
Npt     = 0;
rbounds = [0 1];
name    = 'NA';

% parse varargin
i = 1;
while(i<=length(varargin))
    switch varargin{i}
        case 'pbounds'
            i             = i + 1;
            pbounds       = varargin{i};
            i             = i + 1;
        case 'nMod'
            i             = i + 1;
            nMod          = varargin{i};
            i             = i + 1;
        case 'nRep'
            i             = i + 1;
            nRep          = varargin{i};
            i             = i + 1;
        case 'Npt'
            i             = i + 1;
            Npt          = varargin{i};
            i             = i + 1;
        case 'rbounds'
            i             = i + 1;
            rbounds          = varargin{i};
            i             = i + 1;
        case 'name'
            i             = i + 1;
            name          = varargin{i};
            i             = i + 1;
    end
end

%% Plot CM

% open figure
fh.cm = figure('Name', name); clf;
set(fh.cm,'position',[10 100 800 800],'paperunits','centimeters','Color','w');

% initiate confusion matrix
CM_LL = zeros(nMod);

for count = 1:nRep
    
    % select parameter value for simulation: LB + (UP - LB) .* rand(n)
    b = pbounds(1,1) + (pbounds(2,1)-pbounds(1,1)) .* rand(1,1);
    epsilon = pbounds(1,2) + (pbounds(2,2)-pbounds(1,2)) .* rand(1,1);
    alpha = pbounds(1,3) + (pbounds(2,3)-pbounds(1,3)) .* rand(1,1);
    beta = pbounds(1,4) + (pbounds(2,4)-pbounds(1,4)) .* rand(1,1);
    alpha_c = pbounds(1,5) + (pbounds(2,5)-pbounds(1,5)) .* rand(1,1);
    beta_c = pbounds(1,6) + (pbounds(2,6)-pbounds(1,6)) .* rand(1,1);
    
    % Model 1
    [a, r, pt] = simulate_M1random_v1(T, rbounds, b, mu, Npt);
    [~, ~, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta, alpha_c, beta_c);
    CM_LL(1,:) = CM_LL(1,:) + BEST;
    
    % Model 2
    [a, r, pt] = simulate_M2WSLS_v1(T, rbounds, epsilon, mu, Npt);
    [~, ~, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta, alpha_c, beta_c);
    CM_LL(2,:) = CM_LL(2,:) + BEST;
    
    % Model 3
    [a, r, pt] = simulate_M3RescorlaWagner_v1(T, alpha, beta, mu, rbounds, Npt);
    [~, ~, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta, alpha_c, beta_c);
    CM_LL(3,:) = CM_LL(3,:) + BEST;
    
     % Model 4
    [a, r, pt] = simulate_M4RWCK_v1(T, alpha, beta, alpha_c, beta_c, mu, rbounds, Npt);
    [~, ~, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta, alpha_c, beta_c);
    CM_LL(4,:) = CM_LL(4,:) + BEST;
    
    % Model 5
    [a, r, pt] = simulate_M5CK_v1(T, alpha_c, beta_c, mu, rbounds, Npt);
    [~, ~, BEST] = lik_all_v1(a, r, pt, b, epsilon, alpha, beta, alpha_c, beta_c);
    CM_LL(5,:) = CM_LL(5,:) + BEST;
    
    % open figure
    figure(fh.cm); clf;
    
    % calculate probability
    FM_LL = round(100*CM_LL/sum(CM_LL(1,:)))/100;
    
    % display number/text on scaled image 
    t = imageTextMatrix(FM_LL);
    
    % set font color of values less than 0.3 to white
    set(t(FM_LL'<0.3), 'color', 'w')
    set(t, 'fontsize', 22)
    
    hold on;
    
    % add matrix lines
    addFacetLines(CM_LL);
    
    % add count number as title
    title(['count = ' num2str(count)]);
   
    % set ticks and labels
    set(gca, 'xtick', 1:nMod, 'ytick', 1:nMod, 'fontsize', 18, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')
    
    % update the figure
    drawnow
   
end

% settings
title(sprintf('%s:\np(fit model | simulated model)', name))
set(gcf, 'Position', [311   217   700   600]) 
set(gca, 'fontsize', 20);

end

