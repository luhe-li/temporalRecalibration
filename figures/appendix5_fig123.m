% Appendix 5 - fig 1-3
% individual model prediction of recalibration effect
% 1. causal inference model
% 2. asynchrony-contingent model
% 3. asynchrony-correction model

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

model_slc = [1,3,5];
sub_slc = [1:4, 6:10];
save_fig = 1;

for mm = model_slc
    result_folder = fullfile(dataDir, 'recalibration_models', folders{mm});
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');

    for ss = 1:numel(sub_slc)
        pred{mm, ss} = R{mm, ss}.pred;
        pred_recal(mm, ss, :) = mean(pred{mm, ss}.pss_shift, 2);
    end
end
%% load atheoretical model results

result_folder = fullfile(dataDir, 'atheoretical_models', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

% initialize toj_pss based on the first loaded file
toj_pss = zeros(numel(sub_slc), size(atheo{1}.pred.pss_shift, 2));
for ss = 1:numel(sub_slc)
    toj_pss(ss, :) = atheo{ss}.pred.pss_shift;
end

% calculate group mean and standard error
mean_toj_pss = mean(toj_pss, 1, 'omitnan');
se_toj_pss = std(toj_pss, [], 1, 'omitnan') ./ sqrt(numel(sub_slc));

%% load bootstrap results of atheoretical models

btst_data = load_subject_data(result_folder, sub_slc, 'diag_btst_sub-*');

btst_pss = [];
for ss = 1:numel(sub_slc)
    btst = btst_data{ss};
    for jj = 1:numel(btst.pred)
        btst_pss(jj, :) = btst.pred{jj}.pss_shift;
    end

    for tt = 1:size(btst_pss, 2)
        [lb(ss, tt), ub(ss, tt)] = get95CI(btst_pss(:, tt));
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 6;
titleSZ = 9;
dotSZ = 10;
adaptor_soa = pred{1,1}.adaptor_soa; %ms

for mm = model_slc

    figure; hold on
    set(gcf, 'Position',[1,1, 420, 300]);
    sgtitle([specifications{mm} ' model'],'FontSize',titleSZ)

    for ss = 1:numel(sub_slc)

        subplot(3,3,ss); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
        title(['S' num2str(ss)],'FontSize',titleSZ)

        %% plot atheoretical prediction

        blackL  = plot(adaptor_soa, toj_pss(ss,:),'ko', 'MarkerFaceColor','k','MarkerSize',2);

        for jj = 1:9
            plot([adaptor_soa(jj),adaptor_soa(jj)],[lb(ss, jj), ub(ss, jj)],'k-','LineWidth',lw)
        end

        %% plot model prediction

        redL = plot(adaptor_soa, squeeze(pred_recal(mm, ss,:)), '-o','LineWidth',lw, 'Color','r','MarkerSize',2);

        % look better
        yl = 250;
        ylim([-yl, yl])
        yticks([-yl, 0, yl])
        yticklabels([-yl, 0, yl]./1e3)
        yline(0,'--','LineWidth',1)
        xticks(adaptor_soa)
        xticklabels(adaptor_soa/1e3)
        xtickangle(60)
        xlim([min(adaptor_soa)-100, max(adaptor_soa)+100])
        ax.YTick = [-0.2, -0.1, 0, 0.1, 0.2]*1e3;
        ax.YTickLabel = [-0.2, -0.1, 0, 0.1, 0.2];

        if ss == 4
            ylabel('Recalibration effect (s)','FontSize',titleSZ)
        elseif ss == 8
            xlabel('Adapter SOA (s)','FontSize',titleSZ)
        end

        if ss == numel(sub_slc)
            lgd = legend([blackL, redL],{'Empirical data','Model prediction'});
            ldg.LineWidth = lw;
            lgd.ItemTokenSize = [10,10];
            saveas(gca, fullfile(out_dir, sprintf('M-%s_indiv_recal', folders{mm})),'pdf')
        end

    end

end