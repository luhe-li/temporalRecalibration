% fig 2. Plot recalibration effect from atheoretical model prediction (exponential, shift bias)
% A. Example participant, example session
% B. Group recalibration effect

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load

sub_slc = [1:4,6:10];
exp_sub = 6; % sub7
exp_adp = 2;
result_folder = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');
files = dir(fullfile(result_folder, 'sub-*'));
btst_files = dir(fullfile(result_folder, 'btst_sub-*'));

for ss = 1:numel(sub_slc)

    i_sub = sub_slc(ss);
    i_data = load(fullfile(result_folder, files(i_sub).name));
    pred{ss} = i_data.pred;
    toj_pss(ss,:) = i_data.pred.pss_shift;

    %     % load bootstrap data
    %     btst_data = load(fullfile(result_folder, btst_files(ss).name));
    %
    %     % for each btst trial
    %     for jj = 1:size(btst_fits, 1)
    %         btst_recal(jj, :) = mean(btst_data.pred{jj}.pss_shift, 2);
    %     end
    %
    %     % get confidence interval for each time point
    %     for tt = 1:num_time_points
    %         [toj_lb(tt), toj_ub(tt)] = get68CI(btst_recal(:, tt));
    %     end

end

% calculate group mean
mean_toj_pss   = mean(toj_pss, 1, 'omitnan');
se_toj_pss     = std(toj_pss, [], 1, 'omitnan')./sqrt(numel(sub_slc));

% reorganize data
for ss = 1:numel(sub_slc)
    sub = sub_slc(ss);
    for  ses  = 1:9
        tempD(ses)  = organizeData(sub, ses);
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

adaptor_soa = pred{1,1}.adaptor_soa;

%% A. individual behavior

figure;
set(gcf, 'Position',[0,0,420,130]);

subplot(1,2,1); 
set(gca, 'Position', [0.1, 0.2, 0.45, 0.8]);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
hold on

jitter = (exp_adp < 4) * 10 - (exp_adp >= 4) * 10;
alpha_value           = [repmat(0.4, 1, 3), repmat(0.4, 1, 3)];

% pre
set(gca, 'ColorOrder', cmp1)

% data
scatter(D(exp_sub).data(exp_adp).pre_ms_unique, D(exp_sub).data(exp_adp).pre_pResp, dotSZ,'LineWidth',lw);

% prediction
plot(pred{exp_sub}.test_soa, pred{exp_sub}.pre_pmf{exp_adp},'--', 'LineWidth',lw)

% errorbar

% post
set(gca, 'ColorOrder', cmp2)

% data
scatter(D(exp_sub).data(exp_adp).post_ms_unique, D(exp_sub).data(exp_adp).post_pResp, dotSZ,'filled');

% prediction
plot(pred{exp_sub}.test_soa, pred{exp_sub}.post_pmf{exp_adp},'LineWidth',lw)

% errorbar

% adaptor
xline(adaptor_soa(exp_adp))

% look better
xlim([-550, 550])
x_ticks = [-500, -300:100:300, 500];
y_ticks = [0:0.25:1];
xticks(x_ticks)
xticklabels(strsplit(num2str(x_ticks./1000)))
yticks(y_ticks)
xtickangle(45)
yticklabels(strsplit(num2str(y_ticks)))
xlabel('Test SOA (s)', 'FontName', 'Helvetica', 'FontWeight', 'Light');
ylabel('Probability', 'FontName', 'Helvetica', 'FontWeight', 'Light');

%% B. group behavior

subplot(1,2,2); hold on
set(gca, 'Position', [0.65, 0.2, 0.3, 0.8]);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontWeight', 'Light');
hold on

% group
e = errorbar(adaptor_soa, mean_toj_pss, se_toj_pss, ...
    'o','LineWidth', 1.5);
e.CapSize = 0; e.Color = 'k'; e.MarkerFaceColor = 'k';
e.LineWidth = lw;
e.MarkerSize   = 3.5;

yl = 100;
ylim([-yl, yl])
yticks([-yl, 0, yl])
yticklabels([-yl, 0, yl]./1e3)
yline(0,'--','LineWidth',lw)
ylabel('Recalibration effect (s)','FontName', 'Helvetica', 'FontWeight', 'Light');
xticks(adaptor_soa)
xticklabels(adaptor_soa/1e3)
xtickangle(45)
xlim([min(adaptor_soa)-50, max(adaptor_soa)+50])
xlabel('Adaptor SOA (s)','FontName', 'Helvetica', 'FontWeight', 'Light');
set(gca,'TickDir','out');

flnm = '2TOJ_behavior';
saveas(gca, fullfile(out_dir, flnm),'pdf')
