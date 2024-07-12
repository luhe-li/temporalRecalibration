% fig 4. Simulation of the causal inference model
% A. effect of increasing bias, \beta on recalibration magnitude
% B. effect of increasing sensory noise of one modality, \tau_a or \tau_v,
% on recalibration asymmetry

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

fileName = 'sim_pcc.mat';
if exist(fullfile(out_dir, fileName), 'file') == 2  
    fprintf('File found! Loading %s\n', fileName);
    load(fileName);
else
    fprintf('File does not exist. Performing simulation...\n');

    %% sim pcc
    pccs = linspace(0, 1, 7);
    n_level = numel(pccs);
    parfor i = 1:n_level
        tempModel = currModel;
        pcc = pccs(i);
        pred =  tempModel([beta, tau_a, tau_v, criterion, lambda, pcc, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_a(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fileName)

end

%% simulation B

fileName = 'sim_taus.mat';
if exist(fullfile(out_dir, fileName), 'file') == 2  
    fprintf('File found! Loading %s\n', fileName);
    load(fileName);  
else

    fprintf('File does not exist. Performing simulation...\n');

    %% 1. sim tau_a
    tau_as = [40:5:60];
    n_level = numel(tau_as);
    parfor i = 1:n_level
        tempModel = currModel;
        i_tau_a = tau_as(i);
        pred =  tempModel([beta, i_tau_a, tau_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_a(i,:)       = mean(pred.pss_shift, 2);
    end

    %% 2. sim tau_v
    tau_vs = [40:5:60];
    n_level = numel(tau_vs);
    parfor i = 1:n_level
        tempModel = currModel;
        i_tau_v = tau_vs(i);
        pred =  tempModel([beta, tau_a, i_tau_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
        recal_tau_v(i,:)       = mean(pred.pss_shift, 2);
    end
    save(fileName)
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
    [grad_pcc{c},im_pcc{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end


% color for taus, dark to light
blue = [30, 120, 180; 166, 206, 227]./255;
red = [227, 27, 27; 251, 154, 153]./255;
colorpicked = {blue, red};
depth = n_level;

for c = 1:numel(colorpicked)
    if c ==2
        depth = 5;
    end
    [grad{c},im{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end

%% A. plot p_CC

figure;
set(gcf, 'Position', [0,0,420,150]); hold on

subplot(1,2,1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
set(gca, 'ColorOrder', grad_pcc{1});
% axis equal

% Plot your data with labels
plot(ori_adaptor_soa, recal1,'LineWidth',lw*2);

% Create an array of legend labels corresponding to \tau values
legendLabels = cell(1, numel(pccs));
for i = 1:numel(pccs)
    legendLabels{i} = sprintf('%.2f', pccs(i));
end

% Add the legend with specified labels
leg  = legend(legendLabels, 'Location', 'northwest');
leg.Title.String = 'p_{common}';
leg.ItemTokenSize = [repmat(10,1,7)];

% Add yline (excluding it from the legend)
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

yl = 300;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
ylabel('Recalibration effect (s)')
xlabel('Adaptor SOA (s)')
xticks(ori_adaptor_soa)
xticklabels(ori_adaptor_soa/1e3)
xlim([min(ori_adaptor_soa)-50, max(ori_adaptor_soa)+50])

flnm = 'sim_pcc';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% B1. plot tau_a
figure;
set(gcf, 'Position', [0,0,420,150]); hold on

subplot(1,2,1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
set(gca, 'ColorOrder', grad{1});

% Plot your data with labels
plot(adaptor, recal_tau_a,'LineWidth',lw*2)

% % Create an array of legend labels corresponding to \beta values
% legendLabels = cellstr(num2str(tau_as(:)));
%
% % Add the legend with specified labels
% lgd = legend(legendLabels, 'Location', 'northwest');
% lgd.Title.String = '\tau_{A}';
% ldg.LineWidth = lw;
% lgd.ItemTokenSize = [10,10];

% Add yline (excluding it from the legend)
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

%% B2. plot tau_v

subplot(1,2,2); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out')
set(gca, 'ColorOrder', grad{2});

% Plot your data with labels
plot(adaptor, recal_tau_v,'LineWidth',lw*2)

% % Create an array of legend labels corresponding to \beta values
% legendLabels = cellstr(num2str(tau_vs(:)));
%
% % Add the legend with specified labels
% lgd = legend(legendLabels, 'Location', 'northwest');
% lgd.Title.String = '\tau_{V}';
% ldg.LineWidth = lw;
% lgd.ItemTokenSize = [10,10];

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
xlabel('Adaptor SOA (s)')

flnm = 'sim_tau';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% B1.1-1.2 zoom-in figures of likelihood

figure;
set(gcf, 'Position', [0,0,210,70]); hold on

subplot(1,2,1); hold on
set(gca, 'ColorOrder', grad{1});

soa = 0;
min_x = -200;
max_x = 200;
x_axis = min_x:max_x;
n_level = numel(tau_as);

idx_peak = find(x_axis == (soa + beta));
lx_axis = x_axis(1:idx_peak);
rx_axis = x_axis(idx_peak+1:end);

for i = n_level:-1:1

    i_tau_a = tau_as(i);

    lf = (1/(i_tau_a + tau_v)).* exp(1/i_tau_a * (lx_axis - ((soa + beta))));
    rf = (1/(i_tau_a + tau_v)).* exp(-1/tau_v * (rx_axis - ((soa + beta))));
    likelihood = [lf, rf];
    likelihood = likelihood./sum(likelihood);

    plot(x_axis, likelihood,'-','LineWidth',lw/2);

end

xticks(0)
yticks([])

subplot(1,2,2); hold on
set(gca, 'ColorOrder', grad{2});

soa = 0;
min_x = -200;
max_x = 200;
x_axis = min_x:max_x;
n_level = numel(tau_vs);

idx_peak = find(x_axis == (soa + beta));
lx_axis = x_axis(1:idx_peak);
rx_axis = x_axis(idx_peak+1:end);

for i = n_level:-1:1

    i_tau_v = tau_vs(i);

    lf = (1/(tau_a + i_tau_v)).* exp(1/tau_a * (lx_axis - ((soa + beta))));
    rf = (1/(tau_a + i_tau_v)).* exp(-1/i_tau_v * (rx_axis - ((soa + beta))));
    likelihood = [lf, rf];
    likelihood = likelihood./sum(likelihood);

    plot(x_axis, likelihood,'-','LineWidth',lw/2);

end


xticks(0)
yticks([])

flnm = 'sim_tau_likelihood_demo';
saveas(gca,fullfile(out_dir,flnm),'pdf')

