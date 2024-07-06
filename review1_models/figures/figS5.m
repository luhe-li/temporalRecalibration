% fig S1: model prediction of the causal-inference model with
% modality-specific uncertainty

clear; clc; close all;

%% model info

specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Atheoretical'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','exp_shiftMu'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load recal models

save_fig = 1;
model_slc = 1; 
sub_slc = [1:4,6:10];

recal_folder = fullfile(projectDir, 'recalibration_models_VBMC', folders{model_slc});
files = dir(fullfile(recal_folder, 'sub-*'));

for ss = 1:numel(sub_slc)

    i_sub = sub_slc(ss);
    i_data = load(fullfile(recal_folder, files(ss).name));
    log_model_evi(ss) = i_data.diag.bestELCBO;
    bestP{ss} = i_data.diag.post_mean;
    pred{ss} = i_data.pred;

    % extract summary data for plot
    pred_recal(ss, :) = mean(pred{ss}.pss_shift,2);

end

%% load atheoretical model

athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');
files = dir(fullfile(athe_path, 'sub-*'));

for ss = 1:numel(sub_slc)
    i_sub = sub_slc(ss);
    i_data = load(fullfile(athe_path, files(ss).name));
    toj_pss(ss,:) = i_data.pred.pss_shift;
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

figure; hold on
set(gcf, 'Position',[1,1, 1500, 1000]);
adaptor_soa = pred{1,1}.adaptor_soa; %ms

for ss = 1:numel(sub_slc)

    sub = sub_slc(ss);

%     %% load btst data
%     flnm = sprintf('btst_sub%i_', sub);
%     allFiles = dir(fullfile(athe_path, [flnm '*.mat']));
%     btst_fits = load(allFiles.name);
%     
%     % for each btst trial
%     for jj = 1:size(btst_fits)
%         btst_recal(jj,:) = mean(pred{jj}.pss_shift,2);
%     end
% 
%     % for each time point
%     for tt = 1:size(btst_recal,2)
%         [toj_lb(tt), toj_ub(tt)] = get95CI(btst_recal(:, tt));
%     end

    %% plot
    subplot(3,3,ss); hold on
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    title(['S' num2str(ss)],'FontSize',titleSZ)

    %% plot atheoretical prediction

    l  = plot(adaptor_soa, toj_pss(ss,:),'ko', 'MarkerFaceColor','k');

%     for jj = 1:9
%         plot([adaptor_soa(jj),adaptor_soa(jj)],[-toj_lb(jj), -toj_ub(jj)],'k-','LineWidth',1)
%     end

    %% plot model prediction

    plot(adaptor_soa, pred_recal(ss,:), '-o','LineWidth',lw, 'Color','r')

    % look better
    yl = 250;
    ylim([-yl, yl])
    yticks([-yl, 0, yl])
    yticklabels([-yl, 0, yl]./1e3)
    yline(0,'--','LineWidth',1.5)
    xticks(adaptor_soa)
    xticklabels(adaptor_soa/1e3)
    xtickangle(60)
    xlim([min(adaptor_soa)-100, max(adaptor_soa)+100])
    set(gca,'TickDir','out');

    if ss                              == numel(sub_slc)
        if save_fig
            flnm = 'indiv_recal';
            saveas(gca, fullfile(out_dir, flnm),'pdf')
        end
    end

    if ss == 4
        ylabel('Recalibration effect (s)','FontWeight','bold','FontSize',titleSZ)
    elseif ss == 8
        xlabel('Adapter SOA (s)','FontWeight','bold','FontSize',titleSZ)
    end

end