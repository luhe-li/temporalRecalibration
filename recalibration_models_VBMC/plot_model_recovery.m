clear; close all;

%% select models

specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'};
% specifications = {'Heuristic\newlineModality-specific\newlineuncertainty',...
%     'Heuristic\newlineModality-independent\newlineuncertainty',...
%     'Causal inference\newlineModality-specific\newlineuncertainty',...
%     'Causal inference\newlineModality-indepdent\newlineuncertainty',...
%     'Fixed updated\newlineModality-specific\newlineuncertainty',...
%     'Fixed updated\newlineModality-indepdent\newlineuncertainty'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(cur_dir);
[git_dir, ~] = fileparts(project_dir);
% dataDir = fullfile(git_dir,'temporalRecalibrationData');
dataDir = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Google Drive','My Drive','temporalRecalibrationData');
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(project_dir, 'vbmc')));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load results

results_folder = fullfile(dataDir,'recalibration_models_VBMC','model_recovery_s3');
files = dir(fullfile(results_folder, 'fitM*'));
pattern = 'fitM(\d+)_sample-(\d+)_';

model_slc = [1,2,5,6];
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

        % Find the index of fit_m with the maximum log_model_evidence for the current i_sample
        [~, max_fit_m_index] = max(log_model_evidence(sim_m, :, i_sample));
        % Increment the count for the corresponding fit_m
        CM(sim_m, max_fit_m_index) = CM(sim_m, max_fit_m_index) + 1;

    end
end
CM = CM./num_i_sample;

% use RMSE of the recalibration prediction for CM
RMSE = RMSE(model_slc,model_slc,:);
[num_sim_m, num_fit_m, num_i_sample] = size(RMSE);
CM2 = zeros(num_sim_m, num_fit_m);
for sim_m = 1:num_sim_m

    for i_sample = 1:num_i_sample

        % Find the index of fit_m with the maximum log_model_evidence for the current i_sample
        [~, min_fit_m_index] = min(RMSE(sim_m, :, i_sample));
        % Increment the count for the corresponding fit_m
        CM2(sim_m, min_fit_m_index) = CM2(sim_m, min_fit_m_index) + 1;

    end
end
CM2 = CM2./num_i_sample;

%% plot 

figure;
set(gcf, 'Position',[0,0,420,300]);

imagesc(CM);
colormap('bone')
xticks(1:6)
yticks(1:6)
xticklabels(specifications(model_slc))
yticklabels(specifications(model_slc))
xlabel('Model used for fitting','FontWeight','bold');
ylabel('Data generating model','FontWeight','bold');
title('Percentage of winning based on model evidence','fontsize',12)

[num_rows, num_cols] = size(CM);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM(row, col);
        % Choose text color for better contrast
        textColor = 'w'; % default black
        if val >= 0.3
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end

saveas(gca, fullfile(out_dir, 'CM'),'png');

%% plot2: RMSE of the recalibration phase as the criterion of winning model

figure;
set(gcf, 'Position',[0,0,420,300]);

imagesc(CM2);
colormap('bone')
xticks(1:6)
yticks(1:6)
xticklabels(specifications(model_slc))
yticklabels(specifications(model_slc))
xlabel('Model used for fitting','FontWeight','bold');
ylabel('Data generating model','FontWeight','bold');
title({'Percentage of winning based on';'RMSE of recalibration prediction'},'fontsize',12)

[num_rows, num_cols] = size(CM2);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM2(row, col);
        % Choose text color for better contrast
        textColor = 'w'; % default black
        if val >= 0.3
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end

saveas(gca, fullfile(out_dir, 'CM_RMSE'),'png');

%% compare model predictions
% 
% i_sample = 18;
% sim_m = 1;
% for fit_m = [1,4]
%     figure;
%     set(gcf, 'Position',[0 0 500 400])
%     sgtitle(sprintf('sim-Model %i, fit model%i, logME %.2f', sim_m, fit_m, elbo(sim_m, fit_m, i_sample)))
%     hold on
%     plot(pred{sim_m, fit_m, i_sample}.adaptor_soa, mean(pred{sim_m, fit_m, i_sample}.pss_shift,2),'--or')
%     plot(fake_pred{sim_m, fit_m, i_sample}.adaptor_soa, mean(fake_pred{sim_m, fit_m, i_sample}.pss_shift,2),'k')
%     saveas(gca, fullfile(out_dir, sprintf('recal_sim-Model %i, fit model%i', sim_m, fit_m)),'png');
% 
%     figure;  hold on
%     set(gcf, 'Position',[0 0 1000 400])
%     sgtitle(sprintf('sim-Model %i, fit model%i, logME %.2f', sim_m, fit_m, elbo(sim_m, fit_m, i_sample)))
%     subplot(2,5,1); hold on
%     title('pretest')
%     plot(pred{sim_m, fit_m, i_sample}.test_soa, pred{sim_m, fit_m, i_sample}.pre_pmf','--or')
%     plot(fake_pred{sim_m, fit_m, i_sample}.test_soa, fake_pred{sim_m, fit_m, i_sample}.pre_pmf,'k')
% 
%     for ses = 1:9
%         subplot(2,5,ses+1);
%         hold on
%         plot(pred{sim_m, fit_m, i_sample}.test_soa, squeeze(pred{sim_m, fit_m, i_sample}.post_pmf(ses,:,:))','--or')
%         plot(fake_pred{sim_m, fit_m, i_sample}.test_soa, squeeze((fake_pred{sim_m, fit_m, i_sample}.post_pmf(ses,:,:)))','k')
%         title(['posttest' num2str(ses)])
%     end
%     saveas(gca, fullfile(out_dir, sprintf('toj_sim-Model %i, fit model%i', sim_m, fit_m)),'png');
% end