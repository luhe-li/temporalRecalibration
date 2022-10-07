% a04_p1_model_comparison_v5

% This script runs model comparison for ternary TOJ in the pre/post test,
% obtains AIC from 6 models and the fitted parameters

% fit 6 full models (up to 54 parameters):
% (mu, sigma, criterion, lambda) x (pre, post) x 9 sessions

% It's the same thing as v3 but codes are cleaner

%%
clear all; close all; clc;

% set data path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03expCode')-1);
addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));
addpath('./model_comparison_function_v5')

% subject info
all_sub = 1:10;
n_sub = numel(all_sub);
n_all_sub = all_sub(end);

% model specification and initialization
n_model = 6;

% pre-allocate number of subjects x number of models
[AIC, NLL] = deal(NaN(n_all_sub, n_model));

% pre-allocate struct of fitted parameters
% in M1_to_6 (size: 1 x n_all_sub) each element is a cell (size: num_session x
% num_models), within each cell is a vector of fitted parameters
estP = cell(1, n_all_sub);

% fit models, obtain AIC, NLL, estimated P for each subject
parfor s = all_sub

    % get deltaAIC from 6 full models
    [AIC(s,:) , NLL(s,:), estP{s}] = model_comp_fit_each_session_v5(s);

    disp(s)

end

%% plot delta_AIC

deltaAIC = AIC - min(AIC, [], 2);
figure
heatmap(deltaAIC', 'XLabel','subjects', 'Colormap', flipud(bone),...
    'ColorLimits', [0, 14], 'ColorbarVisible', 'off', 'GridVisible', 'off',...
    'FontSize', 15, 'YLabel', 'model'); colorbar

save('a04_mc_results_v5_n10')
