%Model comparison for ternary TOJ in the pre/post test: v1
% assume 7 free parameters in the full model

clear all; close all; clc;

subjIDs = 3:8;
nsub = length(subjIDs);
deltaAIC = NaN(nsub,5); % number of subjects x number of models
s_estP = cell(1, nsub);
for i = 1:nsub 
    subjID = subjIDs(i);
    [deltaAIC(subjID,:) s_estP{subjID}] = model_comp_fitting(subjID);
end

%% plotting table
figure
heatmap(deltaAIC', 'XLabel','subjects', 'Colormap', flipud(bone),...
    'ColorLimits', [0, 14], 'ColorbarVisible', 'off', 'GridVisible', 'off',...
    'FontSize', 15, 'YLabel', 'model'); colorbar
 
%% plotting curves with the best-fitting parameters
subj = 5;
sess = 4;
bestModel = 1;
bestP = s_estP{subj}{sess, bestModel};
model_comp_plotting(subj, sess, bestModel, bestP)

