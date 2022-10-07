% a04_p1_model_comparison_v5

% This script add individual subject for model comparison for ternary TOJ
% in the pre/post test, obtains AIC from 6 models and the fitted parameters

% fit 6 full models (up to 54 parameters):
% (mu, sigma, criterion, lambda) x (pre, post) x 9 sessions

% It's the same thing as v3 but codes are cleaner

%%
clear all; close all; clc;

% set function path
addpath(genpath('model_comparison_function_v5'))

load('a04_mc_results_v5')
all_sub = 1;

% fit models, obtain AIC, NLL, estimated P for each subject
for s = all_sub

    % get deltaAIC from 6 full models
    [AIC(s,:) , NLL(s,:), estP{s}] = model_comp_fit_each_session_v5(s);

end

%% plot delta_AIC

deltaAIC = AIC - min(AIC, [], 2);
figure
heatmap(deltaAIC', 'XLabel','subjects', 'Colormap', flipud(bone),...
    'ColorLimits', [0, 14], 'ColorbarVisible', 'off', 'GridVisible', 'off',...
    'FontSize', 15, 'YLabel', 'model'); colorbar

save('a04_mc_results_v5')
