% a04_p1_model_comparison_v4

% This script runs model comparison for ternary TOJ in the pre/post test,
% obtains AIC from 8 models and the fitted parameters

% fit 6 full models (up to 54 parameters):
% (mu, sigma, criterion) x (pre, post) x 9 sessions, fix lambda

% fit the reduced model:
% M7 (12 parameters):
% 1 mu_pre, 9 mu_post, 1 sigma, 1 criterion, fix lambda
% M8 (14 parameters):
% 1 mu_pre, 9 mu_post, 1 sigma_pre, 1 sigma_post, 1 criterion_pre, 1 criterion_post, fix lambda


%%
clear all; close all; clc;

% set function path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03ExpCode')-1);
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));
addpath(genpath('model_comparison_function'))

% subject info
all_sub = 3:10;
n_sub = numel(all_sub);
n_all_sub = all_sub(end);

% model specification and initialization
lambda = 0.03;
n_model = 8;

% pre-allocate number of subjects x number of models
[AIC, NLL] = deal(NaN(n_all_sub, n_model));

% pre-allocate struct of fitted parameters 

% In M7_to_8 (size: 1 x n_all_sub) each element is a cell (size: 1 x
% num_models), within each cell is a vector of fitted parameters across
% all sessions;
estP.M7_to_8 = cell(1, n_all_sub);
% in M1_to_6 (size: 1 x n_all_sub) each element is a cell (size: num_session x
% num_models), within each cell is a vector of fitted parameters for each
% session.
estP.M1_to_6 = cell(1, n_all_sub);

% fit models, obtain AIC, NLL, estimated P for each subject
for s = all_sub

    % get AIC and estimated parameters from M7, fitting 9 sessions together
     [AIC(s, 7:8), NLL(s, 7:8), estP.M7_to_8{s}] = model_comp_fit_all_sessions_v4(s, lambda);

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
