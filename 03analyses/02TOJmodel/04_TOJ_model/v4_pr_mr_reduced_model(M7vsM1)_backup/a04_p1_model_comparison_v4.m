% a4_p1_model_comparison_v4

% This script runs model comparison for ternary TOJ in the pre/post test,
% obtains AIC from 8 models and the fitted parameters

% fit 6 full models (up to 54 parameters):
% (mu, sigma, criterion) x (pre, post) x 9 sessions, fix lambda

% fit the reduced model (12 parameters):
% 1 mu_pre, 9 mu_post, 1 sigma, 1 criterion, fix lambda

%%
clear all; close all; clc;

% set function path
addpath(genpath('model_comparison_function'))

% subject info
all_sub = 3:10;
n_sub = numel(all_sub);
n_all_sub = all_sub(end);

% model specification and initialization
lambda = 0.03;
n_model = 7;

% pre-allocate number of subjects x number of models
[AIC, NLL] = deal(NaN(n_all_sub, n_model));

% pre-allocate struct of fitted parameters
% In M7 (size: 1 x n_all_sub) each element is a vector of fitted parameters; 
estP.M7 = cell(1, n_all_sub);
% in M1_to_6 (size: 1 x n_all_sub) each element is a cell (size: num_session x
% num_models), within each cell is a vector of fitted parameters
estP.M1_to_6 = cell(1, n_all_sub);

% fit models, obtain AIC, NLL, estimated P for each subject
for s = all_sub

    % get AIC and estimated parameters from M7, fitting 9 sessions together
    [AIC(s,7), NLL(s,7), estP.M7{s}] = model_comp_fit_all_sessions_v4(s, lambda);

    % get deltaAIC from 6 full models
    [AIC(s,1:6) , NLL(s,1:6), estP.M1_to_6{s}] = model_comp_fit_each_session_v4(s, lambda);

end


%% plot delta_AIC

deltaAIC = AIC - min(AIC, [], 2);
figure
heatmap(deltaAIC', 'XLabel','subjects', 'Colormap', flipud(bone),...
    'ColorLimits', [0, 14], 'ColorbarVisible', 'off', 'GridVisible', 'off',...
    'FontSize', 15, 'YLabel', 'model'); colorbar

save('a04_v4_mc_results')