% fig S1: TOJ model prediction of the causal-inference model with
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
num_ses = 9;

result_folder = fullfile(projectDir, 'recalibration_models_VBMC', folders{model_slc});
files = dir(fullfile(result_folder, 'sub-*'));

for ss = 1:numel(sub_slc)

    i_sub = sub_slc(ss);
    i_data = load(fullfile(result_folder, files(i_sub).name));
    pred{ss} = i_data.pred;

end

%% load atheoretical model

athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');
files = dir(fullfile(athe_path, 'sub-*'));

for ss = 1:numel(sub_slc)
    i_sub = sub_slc(ss);
    D(ss) = load(fullfile(athe_path, files(ss).name)); 
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

tick_y = 0:0.5:1;
tick_x = [-500, 0, 500];
test_soa = pred{1}.test_soa;
adaptor_soa = pred{1}.adaptor_soa;

for  ss  = 1:numel(sub_slc)

    sub              = sub_slc(ss);

    figure; hold on
    set(gcf, 'Position', [1400,0,2500,500]);

    tl = tiledlayout(2,9);
    sgtitle(sprintf('S%i',sub),'FontSize',titleSZ,'FontWeight','bold')

    for adapter = 1:9

        % pre
        nexttile(adapter); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)

        % data
        scatter(D(ss).data(adapter).pre_ms_unique, D(ss).data(adapter).pre_pResp, dotSZ, 'filled');

        % prediction
        plot(pred{ss}.test_soa, pred{ss}.pre_pmf, 'LineWidth',lw)

        % look better
        xlim([-550 550])
        xticks([])
        yticks([])

        if adapter == 1
            ylabel('Pretest probability', 'FontSize',fontSZ,'FontWeight','bold')
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))
        end

        % post
        nexttile(adapter+num_ses); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)

        % data
        scatter(D(ss).data(adapter).post_ms_unique, D(ss).data(adapter).post_pResp, dotSZ, 'filled');

        % prediction
        plot(pred{ss}.test_soa, squeeze(pred{ss}.post_pmf(adapter,:,:)),'LineWidth',lw)
        %             xlabel('Adapter SOA','FontSize',fontSZ,'FontWeight','bold')
        xlabel({'Test SOA',sprintf('Adapter SOA = %.1f s', adaptor_soa(adapter)./1e3)},'FontSize',fontSZ,'FontWeight','bold')

        % look better
        xlim([-550 550])
        xticks(tick_x)
        xticklabels(strsplit(num2str(tick_x./1e3)))
        yticks(tick_y)
        yticklabels(strsplit(num2str(tick_y)))

        if adapter == 1
            ylabel('Posttest probability','FontSize',fontSZ,'FontWeight','bold')
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))
        end

    end

    tl.TileSpacing = 'compact';
    %         tl.XLabel.String = sprintf('S%i: %s',sub,concatenated_str);
    %         tl.XLabel.FontSize = titleSZ;
    %         tl.XLabel.FontWeight = 'bold';
    %         tl.XLabel.String = 'Adapter SOA';
    %         tl.XLabel.FontSize = titleSZ;
    %         tl.XLabel.FontWeight = 'bold';

    if save_fig
        flnm  = sprintf('sub%2d_TOJ_prediction',sub);
        saveas(gca, fullfile(out_dir, flnm),'pdf')
    end


end