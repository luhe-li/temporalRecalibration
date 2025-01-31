% Appendix 1-fig 1
% Model recovery results

clear; close all;

%% select models

specifications = {'Asynchrony-contingent, modality-specific-precision', 'Asynchrony-contingent, modality-independent-precision', 'Causal-inference, modality-specific-precision',  'Causal-inference, modality-independent-precision','Asynchrony-correction, modality-specific-precision', 'Asynchrony-correction, modality-independent-precision','Atheoretical'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','trigger_asym','trigger_sym','exp_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
dataDir = fullfile(projectDir,'fit_results');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load results

results_folder = fullfile(dataDir,'recalibration_models','model_recovery_s3');
files = dir(fullfile(results_folder, 'fitM*'));
pattern = 'fitM(\d+)_sample-(\d+)_';

model_slc = 1:6;
log_model_evidence = nan(6, 6, 100);
for pp = 1:size(files)

    flnm =  files(pp).name;
    r = load(fullfile(results_folder, flnm));
    tokens = regexp(flnm, pattern, 'tokens');
    fit_m = str2double(tokens{1}{1});
    i_sample = str2double(tokens{1}{2});

    % metric fo rmodel comparison
    log_model_evidence(:, fit_m, i_sample) = r.summ.bestELBO;
    RMSE(:, fit_m, i_sample) = r.summ.RMSE;

    % extract model predictions
    fake_pred(:, fit_m, i_sample) = r.fake_pred;
    fake_data(:, fit_m, i_sample) = r.fake_data;
    fit_pred(:,fit_m, i_sample) = r.fit_pred;
    fit_est(:, fit_m,i_sample) = r.fit_est;

end

% use model_evidence for CM
log_model_evidence = log_model_evidence(model_slc,model_slc,:);
[num_sim_m, num_fit_m, num_i_sample] = size(log_model_evidence);
CM = zeros(num_sim_m, num_fit_m);
for sim_m = 1:num_sim_m

    for i_sample = 1:num_i_sample

        [~, max_fit_m_index] = max(log_model_evidence(sim_m, :, i_sample));
        CM(sim_m, max_fit_m_index) = CM(sim_m, max_fit_m_index) + 1;

    end
end
CM = CM./num_i_sample;

%% plot 

figure;
set(gcf, 'Position',[0,0,420,300]);
set(gca, 'FontSize', 4)

imagesc(CM);
colormap('bone')
xticks(1:6)
yticks(1:6)
xticklabels(specifications(model_slc))
yticklabels(specifications(model_slc))
xlabel('Model used for fitting','FontWeight','bold');
ylabel('Data generating model','FontWeight','bold');
title({'Percentage of winning'},'fontsize',12)

[num_rows, num_cols] = size(CM);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM(row, col);
        textColor = 'w'; % default black
        if val >= 0.3
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end

saveas(gca, fullfile(out_dir, 'CM'),'pdf');
