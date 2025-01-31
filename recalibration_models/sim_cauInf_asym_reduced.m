
clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set up model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model.mode       = 'predict';

%% load recalibration model results

model_str = 'cauInf_asym';
currModel = str2func(['nll_' model_str]);
sub_slc = [1:4, 6:10];
save_fig = 1;

result_folder = fullfile(dataDir, 'recalibration_models_VBMC', model_str);
addpath(genpath(fullfile(pwd, model_str)));
R = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)
    p = R{ss}.diag.post_mean;
%     p(end) = 330; % simga_c2
    p(end-1) = 50; % sigma_C1
    bestP{ss} = p;
    pred = currModel(p, model, []);
    pred_recal(ss,:) = mean(pred.pss_shift, 2);

end

%% load atheoretical model results

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

% Initialize toj_pss based on the first loaded file
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

%% load bootstrap results of atheoretical models

btst_data = load_subject_data(result_folder, sub_slc, 'diag_btst_sub-*');

btst_pss = [];
for ss = 1:numel(sub_slc)
    btst = btst_data{ss};
    for jj = 1:numel(btst.pred)
        btst_pss(jj, :) = btst.pred{jj}.pss_shift;
    end

    for tt = 1:size(btst_pss, 2)
        [lb(ss, tt), ub(ss, tt)] = get68CI(btst_pss(:, tt));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;


%% plot individual recal

adaptor_soa = pred.adaptor_soa; %ms

figure; hold on
set(gcf, 'Position',[0, 0, 600, 400]);

sgtitle(['Causal infrence, \sigma_{c1} = 50'],'FontSize',titleSZ)

for ss = 1:numel(sub_slc)

    %% plot
    subplot(3,3,ss); hold on
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    set(gca, 'FontName', 'Helvetica');
    title(['S' num2str(ss)],'FontSize',titleSZ)

    %% plot atheoretical prediction

    blackL  = plot(adaptor_soa, toj_pss(ss,:),'ko', 'MarkerFaceColor','k','MarkerSize',3);

    for jj = 1:9
        plot([adaptor_soa(jj),adaptor_soa(jj)],[lb(ss, jj), ub(ss, jj)],'k-','LineWidth',1)
    end

    %% plot model prediction

    redL = plot(adaptor_soa, squeeze(pred_recal(ss,:)), '-o','LineWidth',lw, 'Color','r','MarkerSize',3);

    % look better
    yl = 250;
    ylim([-yl, yl])
    yline(0,'--','LineWidth',1)
    xticks(adaptor_soa)
    xticklabels(adaptor_soa/1e3)
    xtickangle(60)
    xlim([min(adaptor_soa)-100, max(adaptor_soa)+100])

    ax.YTick = [-0.2, -0.1, 0, 0.1, 0.2]*1e3;
    ax.YTickLabel = [-0.2, -0.1, 0, 0.1, 0.2];

    if ss == 4
        ylabel('Recalibration effect (s)','FontWeight','bold','FontSize',titleSZ)
    elseif ss == 8
        xlabel('Adapter SOA (s)','FontWeight','bold','FontSize',titleSZ)
    end


    if ss == numel(sub_slc)
        legend([blackL, redL],{'Data','Model prediction'})
        saveas(gca, fullfile(out_dir, 'fix_sigma_c1'),'png')
    end

end


