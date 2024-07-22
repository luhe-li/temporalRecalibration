% part 2 of model recovry: organize 100 simulated datesets to 1 file

clear; close all; clc;

%% manage paths

% set cluster option
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

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
    end
end

save(fullfile(outDir, sprintf('sim_data')),'fake_data')