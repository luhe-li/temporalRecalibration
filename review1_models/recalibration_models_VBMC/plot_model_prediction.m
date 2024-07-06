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

model_slc = [1,2];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];
save_fig = 1;

for mm = 1:n_model

    recal_folder = fullfile(pwd, folders{mm});
    files = dir(fullfile(recal_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(recal_folder, files(ss).name));
        DATA(mm, ss) = i_data;
        log_model_evi(mm, ss) = i_data.diag.bestELCBO;
        bestP{mm, ss} = i_data.diag.post_mean;
        pred{mm, ss} = i_data.pred;

        % extract summary data for plot
        pred_recal(mm, ss, :) = mean(pred{mm, ss}.pss_shift,2);

    end
end

%% load atheoretical model

athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');
files = dir(fullfile(athe_path, 'sub-*'));

for ss = 1:numel(sub_slc)
    i_sub = sub_slc(ss);
    i_data = load(fullfile(athe_path, files(ss).name));
    toj_pss(ss,:) = i_data.pred.pss_shift;
end

% calculate group mean
mean_toj_pss   = mean(toj_pss, 1, 'omitnan');
se_toj_pss     = std(toj_pss, [], 1, 'omitnan')./sqrt(numel(sub_slc));


%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;
adaptor_soa = pred{1,1}.adaptor_soa; %ms

%% plot recalibration prediction
figure;
set(gcf, 'Position',[0,0,420,130]);
set(gcf, 'DefaultAxesFontName', 'Helvetica Neue');
set(gcf, 'DefaultAxesFontWeight', 'Light');
set(gcf, 'DefaultTextFontName', 'Helvetica Neue');
set(gcf, 'DefaultTextFontWeight', 'Light');
t = tiledlayout(1,4,'Padding', 'compact', 'TileSpacing', 'compact');

yl = 100;
ytks = {[], [], [], [-yl, 0, yl]};
ytklabels = {[], [], [], [-yl, 0, yl]./1e3};

for mm = 1:n_model

    nexttile; hold on
    %     set(gca, 'Position',[0,0,420,150]);
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    %     axis equal
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    set(gca, 'FontName', 'Helvetica Neue');
    set(gca, 'FontWeight', 'Light');

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

    xlabel(t, 'Adaptor SOA (s)','FontSize',titleSZ);
    ylabel(t,'Recalibration effect (s)','FontSize',titleSZ);

    title(specifications{mm},'FontSize',fontSZ,'FontWeight', 'normal');
    saveas(gca, fullfile(out_dir, sprintf('M-%s_group_recal', folders{mm})),'png')

end

%% plot individual recalibration

for mm = 1:n_model

    figure; hold on
    set(gcf, 'Position',[1,1, 1500, 1000]);
    sgtitle(specifications{mm})

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

        plot(adaptor_soa, squeeze(pred_recal(mm, ss,:)), '-o','LineWidth',lw, 'Color','r')

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

        if ss == 4
            ylabel('Recalibration effect (s)','FontWeight','bold','FontSize',titleSZ)
        elseif ss == 8
            xlabel('Adapter SOA (s)','FontWeight','bold','FontSize',titleSZ)
        end

        if ss                              == numel(sub_slc)
           saveas(gca, fullfile(out_dir, sprintf('M-%s_S%i_recal', folders{mm}, sub)),'png')
        end

    end

end