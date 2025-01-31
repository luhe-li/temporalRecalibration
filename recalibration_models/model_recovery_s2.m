% Step 2 of model recovry: organize 100 simulated datesets to 1 file
% Run locally

clear; close all; clc;

%% manage paths

restoredefaultpath;
[project_dir, ~]= fileparts(pwd);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(git_dir, 'vbmc')));
addpath(genpath(fullfile(project_dir, 'utils')));
outDir = fullfile(project_dir, 'fit_results', 'recalibration_models',mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% load

result_folder = fullfile(project_dir, 'fit_results','recalibration_models', 'model_recovery_s1');
files = dir(fullfile(result_folder, 'sim_data_sample-*'));
pattern = 'sim_data_sample-(\d+)';

for pp = 1:size(files)

    flnm =  files(pp).name;
    r = load(fullfile(result_folder, flnm));
    tokens = regexp(flnm, pattern, 'tokens');
    i_sample = str2double(tokens{1}{1});

    for sim_m = 1:6
        fake_data(sim_m, i_sample).data =  r.fake_data(sim_m,i_sample).data;
        fake_data(sim_m, i_sample).gt_p = r.fake_data(sim_m,i_sample).gt_p;
        fake_data(sim_m, i_sample).pred = r.fake_data(sim_m,i_sample).pred;
        fake_data(sim_m, i_sample).mu_GT = r.fake_data(sim_m,i_sample).mu_GT;
        fake_data(sim_m, i_sample).sd_GT = r.fake_data(sim_m,i_sample).sd_GT;
    end

end

save(fullfile(outDir, sprintf('sim_data')),'fake_data')

% %% optinal: check prediction
% 
% i_sample = 1;
% for sim_m = [1,3,5]
%     figure;
%     sgtitle(['Model' num2str(sim_m)])
%     subplot 121
%     plot(fake_data(sim_m, i_sample).pred.adaptor_soa, mean(fake_data(sim_m, i_sample).pred.pss_shift,2))
%     subplot 122
%     plot(fake_data(sim_m, i_sample).pred.test_soa, fake_data(sim_m, i_sample).pred.pre_pmf')
% end
% 
% %% optinal: check ground-truth
% 
% for sim_m = 1:6
%     clearvars i_gt
%     for i_sample = 1:100
%         i_gt(i_sample,:) = fake_data(sim_m, i_sample).gt_p;
%     end
%     figure; hold on
%     num_p = size(i_gt,2);
%     for pp = 1:num_p
%         subplot(2,3,pp)
%         histogram(i_gt(:,pp),'NumBins',20)
%     end
% end