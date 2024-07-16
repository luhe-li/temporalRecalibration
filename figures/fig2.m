% fig 2. Plot recalibration effect from atheoretical model prediction (exponential, shift bias)
% A. Psychometric function (PMF) of example participant, example session
% B. Group recalibration effect

clear; clc; close all;

%% Manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Load atheoretical model results

sub_slc = [1:4, 6:10];
exp_sub = 7; 
idx_exp_exp = 6;
exp_adp = 2;

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
R = load_subject_data(result_folder, sub_slc, 'sub-*');

for ss = 1:numel(sub_slc)

    pred{ss} = R{ss}.pred;
    toj_pss(ss, :) = R{ss}.pred.pss_shift;

end

% Calculate group mean and standard error
mean_toj_pss = mean(toj_pss, 1, 'omitnan');
se_toj_pss = std(toj_pss, [], 1, 'omitnan') ./ sqrt(numel(sub_slc));

% Reorganize data
D = struct([]);
for ss = 1:numel(sub_slc)
    sub = sub_slc(ss);
    tempD = arrayfun(@(ses) organizeData(sub, ses), 1:9);
    [~, order] = sort([tempD.adaptor_soa]);
    D(ss).data = tempD(order);
end

%% Get confidence interval from bootstrapped data for selected subject and adaptor

% Load bootstrap data
R_btst = load_subject_data(result_folder, exp_sub, 'diag_btst_sub-*');

% Get pmf for each bootstrap trial
for jj = 1:numel(R_btst{1}.pred)
    btst_pre_pmf(jj, :, :) = R_btst{1}.pred{jj}.pre_pmf{exp_adp};
    btst_post_pmf(jj, :, :) = R_btst{1}.pred{jj}.post_pmf{exp_adp};
end

% Get confidence interval for each response type, each time point
pre_lb = zeros(3, size(btst_post_pmf, 3));
pre_ub = zeros(3, size(btst_post_pmf, 3));
post_lb = zeros(3, size(btst_post_pmf, 3));
post_ub = zeros(3, size(btst_post_pmf, 3));

for rr = 1:3
    for tt = 1:size(btst_post_pmf, 3)
        [pre_lb(rr, tt), pre_ub(rr, tt)] = get68CI(squeeze(btst_pre_pmf(:, rr, tt)));
        [post_lb(rr, tt), post_ub(rr, tt)] = get68CI(squeeze(btst_post_pmf(:, rr, tt)));
    end
end

%%  %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure setup
cmp1 = [229, 158, 168; 203, 227, 172; 171, 223, 235] ./ 255;
cmp2 = [216, 49, 91; 175, 213, 128; 88, 193, 238] ./ 255;
lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

adaptor_soa = pred{1}.adaptor_soa;
soa_finer = btst.pred{1}.test_soa;

%% A. Individual behavior
figure;
set(gcf, 'Position', [0, 0, 420, 130]);

subplot(1, 2, 1);
set(gca, 'Position', [0.1, 0.2, 0.45, 0.8], 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica');
hold on;

alpha_value = repmat(0.4, 1, 6);

% Pre
set(gca, 'ColorOrder', cmp1);

% Data
scatter(D(idx_exp_exp).data(exp_adp).pre_ms_unique, D(idx_exp_exp).data(exp_adp).pre_pResp, dotSZ, 'LineWidth', lw);

% Prediction of pmf
plot(pred{idx_exp_exp}.test_soa, pred{idx_exp_exp}.pre_pmf{exp_adp}, '--', 'LineWidth', lw);

% Errorbar
for i = 1:3
    patch([soa_finer, fliplr(soa_finer)], [pre_lb(i, :) fliplr(pre_ub(i, :))], cmp1(i, :), 'EdgeColor', 'none', 'FaceAlpha', alpha_value(i));
end

% Post
set(gca, 'ColorOrder', cmp2);

% Data
scatter(D(idx_exp_exp).data(exp_adp).post_ms_unique, D(idx_exp_exp).data(exp_adp).post_pResp, dotSZ, 'filled');

% Prediction
plot(pred{idx_exp_exp}.test_soa, pred{idx_exp_exp}.post_pmf{exp_adp}, 'LineWidth', lw);

% Errorbar
for i = 1:3
    patch([soa_finer, fliplr(soa_finer)], [post_lb(i, :) fliplr(post_ub(i, :))], cmp2(i, :), 'EdgeColor', 'none', 'FaceAlpha', alpha_value(i));
end

% Adaptor
xline(adaptor_soa(exp_adp));

% Axis settings
xlim([-550, 550]);
x_ticks = [-500, -300:100:300, 500];
y_ticks = 0:0.25:1;
xticks(x_ticks);
xticklabels(x_ticks./1000);
yticks(y_ticks);
xtickangle(60);
yticklabels(y_ticks);
xlabel('Test SOA (s)', 'FontName', 'Helvetica', 'FontWeight', 'Light');
ylabel('Probability', 'FontName', 'Helvetica', 'FontWeight', 'Light');

%% B. Group behavior
subplot(1, 2, 2);
set(gca, 'Position', [0.65, 0.2, 0.3, 0.8], 'LineWidth', lw, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontWeight', 'Light');
hold on;

% Group
e = errorbar(adaptor_soa, mean_toj_pss, se_toj_pss, 'o', 'LineWidth', 1.5);
e.CapSize = 0;
e.Color = 'k';
e.MarkerFaceColor = 'k';
e.LineWidth = lw;
e.MarkerSize = 3.5;

yl = 100;
ylim([-yl, yl]);
yticks([-yl, 0, yl]);
yticklabels([-yl, 0, yl] ./ 1e3);
yline(0, '--', 'LineWidth', lw);
ylabel('Recalibration effect (s)', 'FontName', 'Helvetica', 'FontWeight', 'Light');
xticks(adaptor_soa);
xticklabels(adaptor_soa / 1e3);
xtickangle(60);
xlim([min(adaptor_soa) - 50, max(adaptor_soa) + 50]);
xlabel('Adaptor SOA (s)', 'FontName', 'Helvetica', 'FontWeight', 'Light');
set(gca, 'TickDir', 'out');

% Save figure
flnm = '2TOJ_behavior';
saveas(gcf, fullfile(out_dir, flnm), 'pdf');