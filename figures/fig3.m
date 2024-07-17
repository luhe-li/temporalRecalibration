% fig 3. Plot model predictions of recalibration effect

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
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load recal models

model_slc = [1,2];%[1,2,5];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_models_VBMC', folders{mm});
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');
    
    for ss = 1:numel(sub_slc)
        pred{mm, ss} = R{mm, ss}.pred;
        % Extract summary data for plot
        pred_recal(mm, ss, :) = mean(pred{mm, ss}.pss_shift, 2);
    end
end

%% load atheoretical model

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

toj_pss = zeros(numel(sub_slc), size(atheo{1}.pred.pss_shift, 2));
for ss = 1:numel(sub_slc)
    toj_pss(ss, :) = atheo{ss}.pred.pss_shift;
end

% Calculate group mean and standard error
mean_toj_pss = mean(toj_pss, 1, 'omitnan');
se_toj_pss = std(toj_pss, [], 1, 'omitnan') ./ sqrt(numel(sub_slc));

% Reorganize data
D = struct([]);
for ss = 1:numel(sub_slc)
    for ses = 1:9
        tempD(ses) = organizeData(sub_slc(ss), ses);
    end
    [sorted_adaptor_soa, order] = sort([tempD.adaptor_soa]);
    D(ss).data = tempD(order);
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

%% plot recalibration prediction
figure;
set(gcf, 'Position',[0,0,420,130]);
set(gcf, 'DefaultAxesFontName', 'Helvetica Neue');
set(gcf, 'DefaultAxesFontWeight', 'Light');
set(gcf, 'DefaultTextFontName', 'Helvetica Neue');
set(gcf, 'DefaultTextFontWeight', 'Light');
t = tiledlayout(1,4,'Padding', 'compact', 'TileSpacing', 'compact');

adaptor_soa = pred{1,1}.adaptor_soa; %ms

yl = 100;
ytks = {[], [], [], [-yl, 0, yl]};
ytklabels = {[], [], [], [-yl, 0, yl]./1e3};

for mm = 1:2

    nexttile; hold on
    %     set(gca, 'Position',[0,0,420,150]);
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    %     axis equal
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    set(gca, 'FontName', 'Helvetica');

    e = errorbar(adaptor_soa, mean_toj_pss, se_toj_pss, ...
        'o','LineWidth', 1.5);
    e.CapSize = 0; e.Color = 'k'; e.MarkerFaceColor = 'k';
    e.LineWidth = lw;
    e.MarkerSize   = 2;

    % color parameter for model fitting
    clt = repmat(0.7, 1, 3);

    mean_sim_recal = mean(squeeze(pred_recal(mm,:,:)), 1, 'omitnan');
    se_sim_recal = std(squeeze(pred_recal(mm,:,:)), [], 1, 'omitnan')./sqrt(numel(sub_slc));

    % calculate lower and upper bound using se
    lb_sim_recal  = mean_sim_recal - se_sim_recal;
    ub_sim_recal  = mean_sim_recal + se_sim_recal;

    % plot fitting range
    p = patch([adaptor_soa, fliplr(adaptor_soa)], [lb_sim_recal, fliplr(ub_sim_recal)], ...
        clt, 'FaceAlpha', 0.2, 'EdgeColor','none');

    ylim([-yl, yl])
    yticks(ytks{mm})
    yticklabels(ytklabels{mm})
    yline(0,'--','LineWidth',lw)
    xticks(adaptor_soa)
    xticklabels({'-0.7','-0.3','','','0','','','0.3','0.7'})
    xtickangle(0)
    xlim([min(adaptor_soa)-50, max(adaptor_soa)+50])

    title(specifications{mm},'FontSize',fontSZ,'FontWeight', 'normal');

end

xlabel(t, 'Adaptor SOA (s)','FontSize',titleSZ);
ylabel(t,'Recalibration effect (s)','FontSize',titleSZ);