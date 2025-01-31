% Simulation: when recalibration vanishes depends on bias and uncertainty

clear; clc; close all;
recompute = 1;

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
tau_a = 70;
tau_v = 70;
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
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model.mode       = 'predict';

model_str = 'cauInf_asym';
model.mode       = 'predict';
currModel = str2func(['nll_' model_str]);
addpath(genpath(fullfile(pwd, model_str)));

%% simulate varying modality-independent tau

fileName = 'sim_cauInf_taus.mat';

if ~exist(fullfile(out_dir, fileName), 'file') || recompute

    fprintf('File does not exist. Performing simulation...\n');

    taus = 40:10:100;
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting section

lw = 0.5;
fontsz = 7;
titleFontSz = 10;
adaptor = model.sim_adaptor_soa;

% color for taus, dark to light
blue = [30, 120, 180; 166, 206, 227]./255;
red = [227, 27, 27; 251, 154, 153]./255;
colorpicked = {blue};
depth = n_level;

for c = 1:numel(colorpicked)
    if c ==2
        depth = 5;
    end
    [grad{c},im{c}]= colorGradient(colorpicked{c}(1,:),colorpicked{c}(2,:),depth);
end

%% plot the effect on vanishing point of recalibration

figure;
set(gcf, 'Position', [0,0,420,150]); hold on

subplot(1,2,1); hold on
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
set(gca, 'ColorOrder', grad{1});

plot(adaptor, recal_tau,'LineWidth',lw*2)
yline(0,'--','LineWidth',lw,'HandleVisibility','off')

legendLabels = cellstr(num2str(taus(:)));
lgd = legend(legendLabels, 'Location', 'northwest');
lgd.Title.String = 'Auditory and visual uncertainty';
ldg.LineWidth = lw;
lgd.ItemTokenSize = [10,10];

title('Causal-inference model')
xlabel('Adapter SOA (s)')
yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
ylabel('Recalibration effect (s)')
xticks(adaptor)
xticklabels(adaptor/1e3)
xtickangle(45)
xlim([min(adaptor)-50, max(adaptor)+50])

flnm = 'sim_cauinf_tau';
saveas(gca,fullfile(out_dir,flnm),'pdf')
