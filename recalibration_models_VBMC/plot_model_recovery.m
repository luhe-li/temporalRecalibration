clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(cur_dir);
[git_dir, ~] = fileparts(project_dir);
dataDir = fullfile(git_dir,'temporalRecalibrationData');
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(project_dir, 'vbmc')));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load results

results_folder           = fullfile(dataDir,'recalibration_models_VBMC','model_recovery');
files = dir(fullfile(results_folder, 'fitM*'));
pattern = 'fitM(\d+)_sample-(\d+)_';

fit_m = 1;
for pp = 1:size(files)

    flnm =  files(pp).name;
    r = load(fullfile(results_folder, flnm));
    tokens = regexp(flnm, pattern, 'tokens');
    fit_m = str2double(tokens{1}{1});
    i_sample = str2double(tokens{1}{2});

    log_model_evidence(:, fit_m, i_sample) = r.summ.bestELBO;
    
end

[num_sim_m, num_fit_m, num_i_sample] = size(log_model_evidence);
CM = zeros(num_sim_m, num_fit_m);
for sim_m = 1:num_sim_m
    
    for i_sample = 1:num_i_sample
        % Find the index of fit_m with the maximum log_model_evidence for the current i_sample
        [~, max_fit_m_index] = max(log_model_evidence(sim_m, :, i_sample));
        % Increment the count for the corresponding fit_m
        CM(sim_m, max_fit_m_index) = CM(sim_m, max_fit_m_index) + 1;
    end
end

%% plot

figure; hold on
set(gcf, 'Position',[0,0,170,150]);

imagesc(CM); 
colormap('bone')
xticks(1:2)
yticks(1:2)
xticklabels({'Learning','Percept'})
yticklabels({'Learning','Percept'})
xlabel('Fit Model');
ylabel('Simulated Model');


[num_rows, num_cols] = size(CM);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM(row, col);
        % Choose text color for better contrast
        textColor = 'w'; % default black
        if val > 0.5
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end