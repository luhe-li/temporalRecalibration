% Appendix 2-fig 2
% A. Individual partiicpant's data of recalibration.
% B. Confidence interval of individual asymmetry index

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
dataDir = fullfile(projectDir,'fit_results');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load atheoretical model results

sub_slc = [1:4,6:10];
result_folder = fullfile(dataDir, 'atheoretical_models', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

toj_pss = zeros(numel(sub_slc), size(atheo{1}.pred.pss_shift, 2));
for ss = 1:numel(sub_slc)
    toj_pss(ss, :) = atheo{ss}.pred.pss_shift;
end

% calculate group mean and standard error
mean_toj_pss = mean(toj_pss, 1, 'omitnan');
se_toj_pss = std(toj_pss, [], 1, 'omitnan') ./ sqrt(numel(sub_slc));

%% load bootstrapped results, extract pss_shift and calculate asymmetry index

result_folder = fullfile(dataDir, 'atheoretical_models', 'exp_shiftMu');
atheo_btst = load_subject_data(result_folder, sub_slc, 'diag_btst_sub-*');

for ss = 1:numel(sub_slc)

    % get confidence interval of PSS
    btst = atheo_btst{ss};
    for jj = 1:numel(btst.pred)
        btst_pss(jj, :) = btst.pred{jj}.pss_shift;
        btst_ai(ss, jj)= sum(atheo_btst{ss}.pred{jj}.pss_shift);
    end

    for tt = 1:size(btst_pss, 2)
        [pss_lb(ss, tt), pss_ub(ss, tt)] = get95CI(btst_pss(:, tt));
    end
    [ai_lb(ss), ai_ub(ss)]  = get95CI(btst_ai(ss,:));

end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.5;
fontSZ = 6;
titleSZ = 9;
dotSZ = 10;

%% A. individual recalibration

tick_y = 0:0.5:1;
tick_x = [-500, 0, 500];
test_soa = atheo{1}.model.test_soa;
adaptor_soa = atheo{1}.model.sim_adaptor_soa;

figure;
set(gcf, 'Position', [0,0,420,300]);

for ss = 1:numel(sub_slc)

    subplot(3,3,ss); hold on
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    title(['S' num2str(ss)],'FontSize',titleSZ)

    %% plot atheoretical prediction

    l  = plot(adaptor_soa, toj_pss(ss,:),'ko', 'MarkerFaceColor','k','MarkerSize',2);

    for jj = 1:9
        plot([adaptor_soa(jj),adaptor_soa(jj)],[pss_lb(ss, jj), pss_ub(ss, jj)],'k-','LineWidth',lw)
    end

    % look better
    yl = 250;
    ylim([-yl, yl])
    yticks([-yl, 0, yl])
    yticklabels([-yl, 0, yl]./1e3)
    yline(0,'--','LineWidth',lw)
    xticks(adaptor_soa)
    xticklabels(adaptor_soa/1e3)
    xtickangle(60)
    xlim([min(adaptor_soa)-100, max(adaptor_soa)+100])
    set(gca,'TickDir','out');

    if ss == 4
        ylabel('Recalibration effect (s)','FontSize',titleSZ)
    elseif ss == 8
        xlabel('Adapter SOA (s)','FontSize',titleSZ)
    end

end

flnm = 'A_indiv_recal';
saveas(gca, fullfile(out_dir, flnm),'pdf')

%% B. asymmetry index

figure;
set(gcf, 'Position', [0,0,420,100]);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out');
hold on 

for ss = 1:numel(sub_slc)

    plot([ss,ss],  [ai_lb(ss), ai_ub(ss)],'k','LineWidth',2)
   
end

yline(0,'--')
xlim([0,10])
ylim([-700, 700])
yticks([-700, 0, 700])
yticklabels([-0.7, 0, 0.7])
xticks(1:9)
xticklabels(1:9)
xlabel('Participant','FontSize',titleSZ)
ylabel('Asymmetry index','FontSize',titleSZ)

flnm = 'B_asymmetry_index';
saveas(gca, fullfile(out_dir, flnm),'pdf')
