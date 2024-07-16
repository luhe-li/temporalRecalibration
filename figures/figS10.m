% fig S10. Simulation of the causal inference model
% Effect of increasing bias, \beta on recalibration magnitude

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% free parameters

% initialize with symmetrical parameters
beta = 0;
tau_a = 60;
tau_v = 60;
criterion = 77.23;
lambda = 0.018;
p_common = 0.5;
alpha = 0.0052;
sigma_C1 = 51.9;
sigma_C2 = 261.39;

%% set up model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = numCores; % fit the model multiple times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model.currModelStr = currModelStr; % current model folder

model_str = 'cauInf_asym';
currModel = str2func(['nll_' model_str]);
model.mode       = 'predict';

%% simulation A

fileName = 'sim_beta.mat';
if exist(fullfile(out_dir, fileName), 'file') == 2  
    fprintf('File found! Loading %s\n', fileName);
    load(fullfile(out_dir, fileName));
else
    fprintf('File does not exist. Performing simulation...\n');

    %% sim beta
    betas = [-100:50:100];
    n_level = numel(betas);
    parfor i = 1:n_level
        tempModel = currModel;
        i_beta = betas(i);
        pred =  tempModel([i_beta, tau_a, tau_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_a(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fullfile(out_dir, fileName))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting section

lw = 0.5;
fontsz = 7;
titleFontSz = 10;

% color for pcc, dark to light
purple = [106, 61, 154; 223,209,230]./255;
green = [51, 160, 43; 178, 223, 138]./255;
colorpicked = {purple, green};
depth = 7;
for c = 1:numel(colorpicked)
    if c ==2
        depth = 5;
    end
    [grad{c},im{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end

%% plot beta

subplot(1,2,2); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out')
set(gca, 'ColorOrder', grad{2});

% Plot your data with labels
plot(ori_adaptor_soa, recal2,'LineWidth',lw*2)

% Create an array of legend labels corresponding to \tau values
legendLabels = cell(1, numel(betas));
for i = 1:numel(betas)
    legendLabels{i} = sprintf('%.0f', betas(i));
end

% Add the legend with specified labels
leg  = legend(legendLabels, 'Location', 'northwest');
leg.Title.String = 'beta';
leg.ItemTokenSize = [repmat(10,1,5)];

% Add yline (excluding it from the legend)
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
xlabel('Adaptor SOA (s)')
xticks(ori_adaptor_soa)
xticklabels(ori_adaptor_soa/1e3)
xlim([min(ori_adaptor_soa)-50, max(ori_adaptor_soa)+50])

flnm = 'sim_beta';
saveas(gca,fullfile(out_dir,flnm),'pdf')


