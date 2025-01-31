% Step 1 of model recovry: simulate 100 datesets from 6 models
% - Run locally
% - Model recovery s1-s3 can be easily combined as a single file to run on
%   the cluster

clear; close all; clc;

%% select models

rng('shuffle'); rng('Shuffle');
specifications = {'Asynchrony-contingent, modality-specific-precision', 'Asynchrony-contingent, modality-independent-precision', 'Causal-inference, modality-specific-precision',  'Causal-inference, modality-independent-precision','Asynchrony-correction, modality-specific-precision', 'Asynchrony-correction, modality-independent-precision'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','trigger_asym','trigger_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

num_cores = feature('numcores');
if isempty(gcp('nocreate')); parpool(num_cores-1); end
n_sample = 100; % number of datasets to simulate

%% manage paths

restoredefaultpath;
[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(git_dir, 'vbmc')));
addpath(genpath(fullfile(project_dir, 'utils')));
outDir = fullfile(project_dir, 'fit_results','recalibration_models', mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% organize data
sub_slc = [1:4,6:10];
for sess = 1:9
    data(sess) = organizeData(1, sess);
end

%% define model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 10;
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators

%% sample ground-truth from best parameter estimates

for sim_m = 1:numel(folders)

    sim_str = folders{sim_m};
    sim_func = str2func(['nll_' sim_str]);
    addpath(genpath(fullfile(pwd, sim_str)));

    model.mode = 'initialize';
    Val = sim_func([], model, []);

    % load best estimates
    clearvars bestP
    result_folder = fullfile(project_dir, 'fit_results','recalibration_models', folders{sim_m});
    R = load_subject_data(result_folder, sub_slc, 'sub-*');
    for ss = 1:numel(sub_slc)
        bestP(ss,:) = R{ss}.diag.post_mean;
    end

    % sample 100 ground truth
    mu_GT = mean(bestP, 1);
    sd_GT = std(bestP, [], 1);
    GT_samples = generate_samples(Val, mu_GT, sd_GT, n_sample);
    model.mode = 'predict';

    parfor i_sample = 1:n_sample
        temp_sim_func = sim_func;
        pred =  temp_sim_func(GT_samples, model, []);
        sim_data = simulateData(pred, data);
        fake_data(sim_m, i_sample).data = sim_data;
        fake_data(sim_m, i_sample).gt_p = GT_samples;
        fake_data(sim_m, i_sample).mu_GT = mu_GT;
        fake_data(sim_m, i_sample).sd_GT = sd_GT;
        fake_data(sim_m, i_sample).pred = pred;
    end

end

save(fullfile(outDir, sprintf('sim_data')),'fake_data')

% %% optinal: check predictions
%
% i_sample = 1;
% for sim_m = 1:numel(folders)
%     figure;
%     sgtitle(folders{sim_m})
%     subplot 121
%     plot(fake_data(sim_m, i_sample).pred.adaptor_soa, mean(fake_data(sim_m, i_sample).pred.pss_shift,2))
%     subplot 122
%     plot(fake_data(sim_m, i_sample).pred.test_soa, fake_data(sim_m, i_sample).pred.pre_pmf')
% end
