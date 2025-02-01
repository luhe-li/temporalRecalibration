% Appendix 6 - fig 1-3
% Individual TOJ model prediction of the winning model,
% causal-inference, modality-specific precision model

clear; clc; close all;

%% model info
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

%% load recal models

save_fig = 1;
model_slc = 3;
sub_slc = [1:4,6:10];
num_ses = 9;

result_folder = fullfile(dataDir, 'recalibration_models', folders{model_slc});
R = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)

   pred{ss} = R{ss}.pred;

end

%% load atheoretical model

result_folder = fullfile(dataDir, 'atheoretical_models', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 6;

tick_y = 0:0.5:1;
tick_x = [-500, 0, 500];
test_soa = pred{1}.test_soa;
adaptor_soa = pred{1}.adaptor_soa;

for  ss  = 1:numel(sub_slc)

    figure; hold on
    set(gcf, 'Position', [0,0,800,150]);

    tl = tiledlayout(2,9);
    sgtitle(sprintf('S%i',ss),'FontSize',titleSZ,'FontWeight','bold')

    for adapter = 1:9

        % pre
        nexttile(adapter); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)
        title(sprintf('Adapter SOA\n%.1f s', adaptor_soa(adapter)./1e3),'FontSize',fontSZ-1,'FontWeight','bold')

        % data
        scatter(atheo{ss}.data(adapter).pre_ms_unique, atheo{ss}.data(adapter).pre_pResp, dotSZ, 'filled');

        % prediction
        plot(pred{ss}.test_soa, pred{ss}.pre_pmf, 'LineWidth',lw)

        % look better
        xlim([-550 550])
        xticks([])
        yticks([])

        if adapter == 1
            ylabel({'Pre-test', 'probability'}, 'FontSize',fontSZ,'FontWeight','bold')
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))
        end

        % post
        nexttile(adapter+num_ses); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)
        if adapter == 5
            xlabel('Test SOA','FontSize',fontSZ,'FontWeight','bold')
        end

        % data
        scatter(atheo{ss}.data(adapter).post_ms_unique, atheo{ss}.data(adapter).post_pResp, dotSZ, 'filled');

        % prediction
        plot(pred{ss}.test_soa, squeeze(pred{ss}.post_pmf(adapter,:,:)),'LineWidth',lw)

        % look better
        xlim([-550 550])
        xticks(tick_x)
        xticklabels(strsplit(num2str(tick_x./1e3)))
        yticks(tick_y)
        yticklabels(strsplit(num2str(tick_y)))

        if adapter == 1
            ylabel({'Post-test'; 'probability'},'FontSize',fontSZ,'FontWeight','bold')
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))
        end

    end

    tl.TileSpacing = 'compact';

    if save_fig
        flnm  = sprintf('sub%02d_TOJ_prediction',ss);
        saveas(gca, fullfile(out_dir, flnm),'pdf')
    end

end