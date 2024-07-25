% part 2 of model recovry: organize 100 simulated datesets to 1 file

clear; close all; clc;

%% manage paths

% set cluster option
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end
w
% restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
outDir =  fullfile(projectDir, 'recalibration_models_VBMC', 'model_recovery_s2');
if ~exist(outDir, 'dir'); mkdir(outDir); end
if useCluster == false; projectDir = dataDir; end

%% load

results_folder = fullfile(projectDir, 'recalibration_models_VBMC', 'model_recovery_s1');
files = dir(fullfile(results_folder, 'sim_data_sample-*'));
pattern = 'sim_data_sample-(\d+)';

for pp = 1:size(files)

    flnm =  files(pp).name;
    r = load(fullfile(results_folder, flnm));
    tokens = regexp(flnm, pattern, 'tokens');
    i_sample = str2double(tokens{1}{1});

    for sim_m = 1:6
        fake_data(sim_m, i_sample).data =  r.fake_data(sim_m,i_sample).data;
        fake_data(sim_m, i_sample).gt = r.fake_data(sim_m,i_sample).gt;
        fake_data(sim_m, i_sample).pred = r.fake_data(sim_m,i_sample).pred;
    end

end

save(fullfile(outDir, sprintf('sim_data')),'fake_data')

%% optinal: check prediction

i_sample = 1;
for sim_m = [1,3,5]
    figure;
    sgtitle(['Model' num2str(sim_m)])
    subplot 121
    plot(fake_data(sim_m, i_sample).pred.adaptor_soa, mean(fake_data(sim_m, i_sample).pred.pss_shift,2))
    subplot 122
    plot(fake_data(sim_m, i_sample).pred.test_soa, fake_data(sim_m, i_sample).pred.pre_pmf')
end

%% optinal: check ground-truth
for sim_m = 1:6
    clearvars i_gt
    for i_sample = 1:100
        i_gt(i_sample,:) = fake_data(sim_m, i_sample).gt;
    end
    figure; hold on
    num_p = size(i_gt,2);
    for pp = 1:num_p
        subplot(2,3,pp)
        histogram(i_gt(:,pp),'NumBins',20)
    end
end