% fig 5. Simulation of the causal inference model
% A. Effect of increasing p_common on recalibration magnitude
% B. Effect of increasing sensory noise on recalibration magnitude at large SOA
% C. Effect of increasing sensory noise of one modality, \tau_a or \tau_v,
%    on recalibration asymmetry

clear; clc; close all;
recompute = 0;
model_str = 'cauInf_asym';

%% manage paths

restoredefaultpath;
[projectDir, ~]= fileparts(pwd);
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'recalibration_models', model_str)));
out_dir = fullfile(pwd, mfilename);
addpath(genpath(out_dir))
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% free parameters

% initialize with symmetrical parameters
beta = 0;
tau_a = 80;
tau_v = 80;
criterion = 77.23;
lambda = 0.018;
p_common = 0.7;
alpha = 0.0052;
sigma_C1 = 51.9;
sigma_C2 = 261.39;

%% set up model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 3; % fit the model multiple times, each with a different initialization
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators

currModel = str2func(['nll_' model_str]);
model.mode = 'predict';

%% simulation A, varying probability of a common cause

fileName = 'sim_pcc.mat';
if ~exist(fullfile(out_dir, fileName), 'file') || recompute

    fprintf('File does not exist. Performing simulation...\n');

    pccs = linspace(0, 1, 7);
    n_level = numel(pccs);
    parfor i = 1:n_level
        tempModel = currModel;
        pcc = pccs(i);
        pred =  tempModel([beta, tau_a, tau_v, criterion, lambda, pcc, alpha, sigma_C1, sigma_C2], model, []);
        recal_pcc(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fullfile(out_dir, fileName))

else
    fprintf('File found! Loading %s\n', fileName);
    load(fullfile(out_dir, fileName));
end

%% simulation B, varying uncertainty of both modalities

fileName = 'sim_taus.mat';
if ~exist(fullfile(out_dir, fileName), 'file') || recompute

    fprintf('File does not exist. Performing simulation...\n');

    taus = 40:10:90;
    n_level = numel(taus);
    parfor i = 1:n_level
        tempModel = currModel;
        tau = taus(i);
        pred =  tempModel([beta, tau, tau, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fullfile(out_dir, fileName))

else
    fprintf('File found! Loading %s\n', fileName);
    load(fullfile(out_dir, fileName));
end

%% simulation C, varying uncertainty of one modality

fileName = 'sim_one_tau.mat';
    
if ~exist(fullfile(out_dir, fileName), 'file') || recompute

    fprintf('File does not exist. Performing simulation...\n');

    %% 1. sim tau_a
    tau_as = 40:10:80;
    n_level = numel(tau_as);
    parfor i = 1:n_level
        tempModel = currModel;
        i_tau_a = tau_as(i);
        pred =  tempModel([beta, i_tau_a, tau_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_a(i,:)       = mean(pred.pss_shift, 2);
    end

    %% 2. sim tau_v
    tau_vs = 40:10:80;
    n_level = numel(tau_vs);
    parfor i = 1:n_level
        tempModel = currModel;
        i_tau_v = tau_vs(i);
        pred =  tempModel([beta, tau_a, i_tau_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_v(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fullfile(out_dir, fileName))

else
    fprintf('File found! Loading %s\n', fileName);
    load(fullfile(out_dir, fileName)); 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting section

lw = 0.5;
fontsz = 7;
titleFontSz = 10;
adaptor = model.sim_adaptor_soa;

% color for pcc, dark to light
purple = [106, 61, 154; 230, 230, 230]./255;
green = [51, 160, 43; 230, 230, 230]./255;
colorpicked = {purple, green};
depth = 7;
for c = 1:numel(colorpicked)
    [grad_lin{c},im_lin{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end

% color for taus, dark to light
blue = [30, 120, 180; 230, 230, 230]./255;
red = [227, 27, 27; 230, 230, 230]./255;
colorpicked = {blue, red};
depth = n_level+1;

for c = 1:numel(colorpicked)
    [grad{c},im{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end

%% A. plot p_CC

figure;
set(gcf, 'Position', [0,0,420,150]); hold on

subplot(1,2,1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
set(gca, 'ColorOrder', grad_lin{1});

plot(adaptor, recal_pcc,'LineWidth',lw*2);

% Create an array of legend labels corresponding to \pcc values
legendLabels = cell(1, numel(pccs));
for i = 1:numel(pccs)
    legendLabels{i} = sprintf('%.2f', pccs(i));
end
leg  = legend(legendLabels, 'Location', 'northwest');
leg.Title.String = 'p_{common}';
leg.ItemTokenSize = [repmat(10,1,7)];

% Add yline (excluding it from the legend)
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
ylabel('Recalibration effect (s)')
xlabel('Adapter SOA (s)')
xticks(adaptor)
xticklabels(adaptor/1e3)
xlim([min(adaptor)-50, max(adaptor)+50])

%% B. plot tau

subplot(1,2,2); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out')
set(gca, 'ColorOrder', grad_lin{2});

plot(adaptor, recal_tau,'LineWidth',lw*2);

% Create an array of legend labels corresponding to \taus values
legendLabels = cell(1, numel(taus));
for i = 1:numel(taus)
    legendLabels{i} = sprintf('%.2f', taus(i)/1000);
end
leg  = legend(legendLabels, 'Location', 'northwest');
leg.Title.String = 'Auditory and visual uncertainty';
leg.ItemTokenSize = [repmat(10,1,7)];
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
xlabel('Adapter SOA (s)')
xticks(adaptor)
xticklabels(adaptor/1e3)
xlim([min(adaptor)-50, max(adaptor)+50])

flnm = 'AB_nonlinearity';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% C1. plot tau_a
figure;
set(gcf, 'Position', [0,0,420,150]); hold on

subplot(1,2,1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
set(gca, 'ColorOrder', grad{1});

plot(adaptor, recal_tau_a,'LineWidth',lw*2)
legendLabels = cellstr(num2str(tau_as(:)/1e3));
lgd = legend(legendLabels, 'Location', 'northwest');
lgd.Title.String = '\tau_{A}';
ldg.LineWidth = lw;
lgd.ItemTokenSize = [10,10];

yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
ylabel('Recalibration effect (s)')
xticks(adaptor)
xticklabels(adaptor/1e3)
xtickangle(45)
xlim([min(adaptor)-50, max(adaptor)+50])

%% C2. plot tau_v

subplot(1,2,2); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out')
set(gca, 'ColorOrder', grad{2});

% Plot your data with labels
plot(adaptor, recal_tau_v,'LineWidth',lw*2)

% Create an array of legend labels corresponding to \beta values
legendLabels = cellstr(num2str(tau_vs(:)/1e3));

% Add the legend with specified labels
lgd = legend(legendLabels, 'Location', 'northwest');
lgd.Title.String = '\tau_{V}';
ldg.LineWidth = lw;
lgd.ItemTokenSize = [10,10];

% Add yline (excluding it from the legend)
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([])
xticks(adaptor)
xticklabels(adaptor/1e3)
xtickangle(45)
xlim([min(adaptor)-50, max(adaptor)+50])
xlabel('Adapter SOA (s)')

flnm = 'C_asymmetry';
saveas(gca,fullfile(out_dir,flnm),'pdf')
