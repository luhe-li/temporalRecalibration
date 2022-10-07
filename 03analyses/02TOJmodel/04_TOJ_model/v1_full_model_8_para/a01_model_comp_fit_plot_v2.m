%Model comparison for ternary TOJ in the pre/post test

clear all; close all; clc;

subjIDs = 3:8;
nsub = length(subjIDs);
deltaAIC = NaN(nsub,5); % number of subjects x number of models
s_estP = cell(1, nsub);
for i = 1:nsub
    subjID = subjIDs(i);
    [deltaAIC(subjID,:) s_estP{subjID}] = model_comp_fitting_v2(subjID);
end

%% plotting table
figure
heatmap(deltaAIC', 'XLabel','subjects', 'Colormap', flipud(bone),...
    'ColorLimits', [0, 14], 'ColorbarVisible', 'off', 'GridVisible', 'off',...
    'FontSize', 15, 'YLabel', 'model'); colorbar

%% plotting curves with the best-fitting parameters
nsess = 9;
for i = 1:nsub
    subjID = subjIDs(i);
    for sess = 1:nsess
        bestModel = 1;
        bestP = s_estP{subj}{sess, bestModel};
        model_comp_plotting_v2(subj, sess, bestModel, bestP)
    end
end

